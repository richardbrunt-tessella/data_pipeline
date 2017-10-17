import logging

from tqdm import tqdm

from mrtarget.Settings import Config
from mrtarget.common import Actions
from mrtarget.common import TqdmToLogger
from mrtarget.common.DataStructure import JSONSerializable, denormDict
from mrtarget.common.ElasticsearchLoader import Loader
from mrtarget.common.ElasticsearchQuery import ESQuery
from mrtarget.common.LookupHelpers import LookUpDataRetriever, LookUpDataType
from mrtarget.common.Redis import RedisQueue, RedisQueueWorkerProcess, RedisQueueStatusReporter
from mrtarget.common.Scoring import ScoringMethods, HarmonicSumScorer
from mrtarget.modules.EFO import EFO
from mrtarget.modules.EvidenceString import Evidence, ExtendedInfoGene, ExtendedInfoEFO
from mrtarget.modules.GeneData import Gene
from mrtarget.modules.HPA import HPAExpression, hpa2tissues

logger = logging.getLogger(__name__)
tqdm_out = TqdmToLogger(logger,level=logging.INFO)


class AssociationActions(Actions):
    EXTRACT = 'extract'
    PROCESS = 'process'
    UPLOAD = 'upload'

class AssociationMap(JSONSerializable):
    '''
    stores al the details for "target to disease" and "disease to targets" links
    as  overall or broken down by datasource and datatype
    '''

    def __init__(self,
                 id,
                 label=None,
                 **kwargs):
        self.id = id
        self.label = label
        self._generate_data_container()
        self.__dict__.update(**kwargs)

    def _generate_data_container(self):
        self.overall = set()
        self.datatype = {}
        self.datasource = {}

    def add_datasource(self, ds, entity):
        if ds not in self.datasource:
            self.datasource[ds]=set()
        self.datasource[ds].add(entity)
        dt =Config.DATASOURCE_TO_DATATYPE_MAPPING[ds]
        if dt not in self.datatype:
            self.datatype[dt]=set()
        self.datatype[dt].add(entity)
        self.overall.add(entity)

    def encode_all_vectors(self, id_order):
        self.overall_vector = self.encode_vector(self.overall, id_order)
        self.datatype_vector = {}
        for dt in self.datatype:
            self.datatype_vector[dt] = self.encode_vector(self.datatype[dt], id_order)
        self.datasource_vector = {}
        for ds in self.datasource:
            self.datasource_vector[ds] = self.encode_vector(self.datasource[ds], id_order)

    @staticmethod
    def encode_vector(v, labels):
        return ' '.join([i if i in v else 'n'+i for i in labels])



class AssociationScore(JSONSerializable):

    def __init__(self):
        self.init_scores()

    def init_scores(self):
        """init scores to 0.0 and map to 2-maps init with 0.0"""
        self.overall = 0.0
        (self.datasources, self.datatypes) = denormDict(
            Config.DATASOURCE_TO_DATATYPE_MAPPING, (0.0, 0.0))


class Association(JSONSerializable):

    def __init__(self, target, disease, is_direct):
        self.target = {'id': target}
        self.disease = {'id': disease}
        self.is_direct = is_direct
        self.set_id()

        for method_key, method in ScoringMethods.__dict__.items():
            if not method_key.startswith('_'):
                self.set_scoring_method(method, AssociationScore())

        self.evidence_count = dict(total=0.0,
                                   datatypes={},
                                   datasources={})

        (self.evidence_count['datasources'],
         self.evidence_count['datatypes']) = denormDict(
            Config.DATASOURCE_TO_DATATYPE_MAPPING, (0.0, 0.0))

        self.private = {}
        self.private['facets'] = dict(datatype=[],
                                      datasource=[],
                                      free_text_search=[],
                                      expression_tissues=[])

    def get_scoring_method(self, method):
        if method not in ScoringMethods.__dict__.values():
            raise AttributeError("method need to be a valid ScoringMethods")
        return self.__dict__[method]

    def set_scoring_method(self, method, score):
        if method not in ScoringMethods.__dict__.values():
            raise AttributeError("method need to be a valid ScoringMethods")
        if not isinstance(score, AssociationScore):
            raise AttributeError("score need to be an instance"
                                 "of AssociationScore")
        self.__dict__[method] = score

    def set_id(self):
        self.id = '%s-%s' % (self.target['id'], self.disease['id'])

    def set_target_data(self, gene):
        """get generic gene info"""
        pathway_data = dict(pathway_type_code=[],
                            pathway_code=[])

        GO_terms = dict(biological_process=[],
                        cellular_component=[],
                        molecular_function=[],
                        )

        target_class = dict(level1=[],
                            level2=[])

        uniprot_keywords = []
        #TODO: handle domains
        genes_info=ExtendedInfoGene(gene)
        '''collect data to use for free text search'''

        for el in ['geneid', 'name', 'symbol']:
            self.private['facets']['free_text_search'].append(
                genes_info.data[el])

        if 'facets' in gene._private and 'reactome' in gene._private['facets']:
            pathway_data['pathway_type_code'].extend(gene._private['facets']['reactome']['pathway_type_code'])
            pathway_data['pathway_code'].extend(gene._private['facets']['reactome']['pathway_code'])
        if gene.go:
            for item in gene.go:
                go_code, data = item['id'], item['value']
                try:
                    category,term = data['term'][0], data['term'][2:]
                    if category =='P':
                        GO_terms['biological_process'].append(dict(code=go_code,
                                                                   term=term))
                    elif category =='F':
                        GO_terms['molecular_function'].append(dict(code=go_code,
                                                                   term=term))
                    elif category =='C':
                        GO_terms['cellular_component'].append(dict(code=go_code,
                                                                   term=term))
                except:
                    pass

        if gene.uniprot_keywords:
            uniprot_keywords = gene.uniprot_keywords

        if genes_info:
            self.target[ExtendedInfoGene.root] = genes_info.data

        if pathway_data['pathway_code']:
            pathway_data['pathway_type_code']=list(set(pathway_data['pathway_type_code']))
            pathway_data['pathway_code']=list(set(pathway_data['pathway_code']))
        if 'chembl' in gene.protein_classification and gene.protein_classification['chembl']:
            target_class['level1'].append([i['l1'] for i in gene.protein_classification['chembl'] if 'l1' in i])
            target_class['level2'].append([i['l2'] for i in gene.protein_classification['chembl'] if 'l2' in i])

        '''Add private objects used just for indexing'''

        if pathway_data['pathway_code']:
            self.private['facets']['reactome']= pathway_data
        if uniprot_keywords:
            self.private['facets']['uniprot_keywords'] = uniprot_keywords
        if GO_terms['biological_process'] or \
            GO_terms['molecular_function'] or \
            GO_terms['cellular_component'] :
            self.private['facets']['go'] = GO_terms
        if target_class['level1']:
            self.private['facets']['target_class'] = target_class

    def set_hpa_data(self, hpa):
        '''set a compat hpa expression data into the score object'''
        filteredHPA = hpa2tissues(hpa)
        if filteredHPA is not None and len(filteredHPA) > 0:
            self.private['facets']['expression_tissues'] = filteredHPA

    def set_disease_data(self, efo):
        """get generic efo info"""
        efo_info=ExtendedInfoEFO(efo)
        '''collect data to use for free text search'''
        self.private['facets']['free_text_search'].append(efo_info.data['efo_id'])
        self.private['facets']['free_text_search'].append(efo_info.data['label'])
        self.private['facets']['free_text_search'].extend(efo_info.data['therapeutic_area']['labels'])

        if efo_info:
            self.disease[ExtendedInfoEFO.root] = efo_info.data

    def set_available_datasource(self, ds):
        if ds not in self.private['facets']['datasource']:
            self.private['facets']['datasource'].append(ds)
            self.private['facets']['free_text_search'].append(ds)
    def set_available_datatype(self, dt):
        if dt not in self.private['facets']['datatype']:
            self.private['facets']['datatype'].append(dt)
            self.private['facets']['free_text_search'].append(dt)

    def __bool__(self):
        return self.get_scoring_method(ScoringMethods.HARMONIC_SUM).overall != 0
    def __nonzero__(self):
        return self.__bool__()


class EvidenceScore():
    def __init__(self,
                 evidence_string = None,
                 score= None,
                 datatype = None,
                 datasource = None,
                 is_direct = None):
        if evidence_string is not None:
            e = Evidence(evidence_string).evidence
            self.score = e['scores']['association_score']
            self.datatype = e['type']
            self.datasource = e['sourceID']
        if score is not None:
            self.score = score
        if datatype is not None:
            self.datatype = datatype
        if datasource is not None:
            self.datasource = datasource
        self.is_direct = is_direct


class Scorer():

    def __init__(self):
        pass

    def score(self,target, disease, evidence_scores, is_direct, method  = None):

        ass = Association(target, disease, is_direct)

        # set evidence counts
        for e in evidence_scores:
            # make sure datatype is constrained
            if all([e.datatype in ass.evidence_count['datatypes'],
                    e.datasource in ass.evidence_count['datasources']]):
                ass.evidence_count['total']+=1
                ass.evidence_count['datatypes'][e.datatype]+=1
                ass.evidence_count['datasources'][e.datasource]+=1

                # set facet data
                ass.set_available_datatype(e.datatype)
                ass.set_available_datasource(e.datasource)

        # compute scores
        if (method == ScoringMethods.HARMONIC_SUM) or (method is None):
            self._harmonic_sum(evidence_scores, ass, scale_factor=2)
        if (method == ScoringMethods.SUM) or (method is None):
            self._sum(evidence_scores, ass)
        if (method == ScoringMethods.MAX) or (method is None):
            self._max(evidence_scores, ass)

        return ass

    def _harmonic_sum(self, evidence_scores, ass, max_entries = 100, scale_factor = 1):
        har_sum_score = ass.get_scoring_method(ScoringMethods.HARMONIC_SUM)
        datasource_scorers = {}
        for e in evidence_scores:
            if e.datasource not in datasource_scorers:
                datasource_scorers[e.datasource]= HarmonicSumScorer(buffer=max_entries)
            datasource_scorers[e.datasource].add(e.score)
        '''compute datasource scores'''
        overall_scorer = HarmonicSumScorer(buffer=max_entries)
        for datasource in datasource_scorers:
            har_sum_score.datasources[datasource]=datasource_scorers[datasource].score(scale_factor=scale_factor, cap=1)
            overall_scorer.add(har_sum_score.datasources[datasource])
        '''compute datatype scores'''
        datatypes_scorers = dict()
        for ds in har_sum_score.datasources:
            dt = Config.DATASOURCE_TO_DATATYPE_MAPPING[ds]
            if dt not in datatypes_scorers:
                datatypes_scorers[dt]= HarmonicSumScorer(buffer=max_entries)
            datatypes_scorers[dt].add(har_sum_score.datasources[ds])
        for datatype in datatypes_scorers:
            har_sum_score.datatypes[datatype]=datatypes_scorers[datatype].score(scale_factor=scale_factor)
        '''compute overall scores'''
        har_sum_score.overall = overall_scorer.score(scale_factor=scale_factor)

        return ass

    def _sum(self, evidence_scores, ass):
        sum_score = ass.get_scoring_method(ScoringMethods.SUM)
        for e in evidence_scores:
            sum_score.overall+=e.score
            sum_score.datatypes[e.datatype]+=e.score
            sum_score.datasources[e.datasource]+=e.score

        return

    def _max(self, evidence_score, ass):
        max_score = ass.get_scoring_method(ScoringMethods.MAX)
        for e in evidence_score:
            if e.score > max_score.datasources[e.datasource]:
                max_score.datasources[e.datatype] = e.score
                if e.score > max_score.datatypes[e.datatype]:
                    max_score.datatypes[e.datatype]=e.score
                    if e.score > max_score.overall:
                        max_score.overall=e.score

        return




class TargetDiseaseEvidenceProducer(RedisQueueWorkerProcess):

    def __init__(self,
                 target_q,
                 r_path,
                 target_disease_pair_q,
                 ):
        super(TargetDiseaseEvidenceProducer, self).__init__(queue_in=target_q,
                                                            redis_path=r_path,
                                                            queue_out=target_disease_pair_q)

    def process(self, data):
        target = data
        target_association_map = AssociationMap(target)
        available_evidence = self.es_query.count_evidence_for_target(target)
        if available_evidence:
            self.init_data_cache()
            evidence_iterator = self.es_query.get_evidence_for_target_simple(target, available_evidence)
            # for evidence in tqdm(evidence_iterator,
            #                    desc='fetching evidence for target %s'%target,
            #                    unit=' evidence',
            #                    unit_scale=True,
            #                    total=available_evidence):
            for evidence in evidence_iterator:
                for efo in evidence['private']['efo_codes']:
                    key = (evidence['target']['id'], efo)
                    if key not in self.data_cache:
                        self.data_cache[key] = []
                    row = EvidenceScore(
                        score=evidence['scores']['association_score'] * Config.SCORING_WEIGHTS[evidence['sourceID']],
                        datatype=Config.DATASOURCE_TO_DATATYPE_MAPPING[evidence['sourceID']],
                        datasource=evidence['sourceID'],
                        is_direct=efo == evidence['disease']['id'])
                    self.data_cache[key].append(row)
                target_association_map.add_datasource([evidence['sourceID']],evidence['disease']['id'])

            self.produce_pairs()

    def init_data_cache(self,):
        try:
            del self.data_cache
        except: pass
        self.data_cache = dict()

    def produce_pairs(self):
        for key,evidence in self.data_cache.items():
            is_direct = False
            for e in evidence:
                if e.is_direct:
                    is_direct = True
                    break
            self.put_into_queue_out((key[0],key[1], evidence, is_direct))
        self.init_data_cache()

    def init(self):
        super(TargetDiseaseEvidenceProducer, self).init()
        self.es_query = ESQuery()

    def close(self):
        super(TargetDiseaseEvidenceProducer, self).close()



class TargetMapStorer(RedisQueueWorkerProcess):
    '''
    this storer takes all the targets to disease maps
    store them in elasticsearch
    collect the diseases the targets are linked to 
    put the diseases in a queue for the inverse map to be processed
    '''

    def __init__(self,
                 target_map_q,
                 r_path,
                 disease_q,
                 lookup_data,
                 chunk_size = 1e4,
                 dry_run = False,
                 ):
        super(TargetMapStorer, self).__init__(queue_in=target_map_q,
                                            redis_path=r_path,
                                            queue_out=disease_q)
        self.lookup_data = lookup_data
        self.chunk_size = chunk_size
        self.dry_run = dry_run
        self.loader = None
        self.diseases = set()

    def init(self):
        super(TargetMapStorer, self).init()
        self.lookup_data.set_r_server(self.r_server)
        self.loader = Loader(chunk_size=self.chunk_size,
                             dry_run=self.dry_run)

    def close(self):
        for disease in self.diseases:
            self.put_into_queue_out(disease)
        self.queue_out.set_submission_finished()
        super(TargetMapStorer, self).close()
        self.loader.close()

    def process(self, data):
        as_map = data
        as_map.encode_all_vectors()
        self.diseases|=as_map.overall
        self.loader.put(Config.ELASTICSEARCH_MAP_ASSOCIATION_INDEX_NAME,
                        Config.ELASTICSEARCH_MAP_ASSOCIATION_TARGET_DOC_NAME,
                        as_map.id,
                        as_map
                        )

class DiseaseMapStorer(RedisQueueWorkerProcess):
    '''
    this storer takes all the diseases
    builds disease to target maps
    stores them in elasticsearch
    '''

    def __init__(self,
                 disease_q,
                 r_path,
                 queue_out,
                 lookup_data,
                 chunk_size = 1e4,
                 dry_run = False,
                 ):
        super(DiseaseMapStorer, self).__init__(queue_in=disease_q,
                                            redis_path=r_path,
                                            queue_out=queue_out)
        self.lookup_data = lookup_data
        self.chunk_size = chunk_size
        self.dry_run = dry_run
        self.loader = None

    def init(self):
        super(DiseaseMapStorer, self).init()
        self.lookup_data.set_r_server(self.r_server)
        self.loader = Loader(chunk_size=self.chunk_size,
                             dry_run=self.dry_run)

    def close(self):

        super(DiseaseMapStorer, self).close()
        self.loader.close()

    def process(self, data):
        disease_id = data
        as_map = AssociationMap(disease_id)
        evidence_iterator = self.es_query.get_evidence_for_disease_simple(disease_id)

        for evidence in evidence_iterator:
            as_map.add_datasource([evidence['sourceID']], evidence['target']['id'])
        as_map.encode_all_vectors()
        self.diseases|=as_map.overall
        self.loader.put(Config.ELASTICSEARCH_MAP_ASSOCIATION_INDEX_NAME,
                        Config.ELASTICSEARCH_MAP_ASSOCIATION_DISEASE_DOC_NAME,
                        as_map.id,
                        as_map
                        )

class ScoreProducer(RedisQueueWorkerProcess):

    def __init__(self,
                 evidence_data_q,
                 r_path,
                 score_q,
                 lookup_data,
                 chunk_size = 1e4,
                 dry_run = False,
                 es = None
                 ):
        super(ScoreProducer, self).__init__(queue_in=evidence_data_q,
                                            redis_path=r_path,
                                            queue_out=score_q)
        self.evidence_data_q = evidence_data_q
        self.score_q = score_q
        self.lookup_data = lookup_data
        self.chunk_size = chunk_size
        self.dry_run = dry_run
        self.es = None
        self.loader = None

    def init(self):
        super(ScoreProducer, self).init()
        self.scorer = Scorer()
        self.lookup_data.set_r_server(self.r_server)
        self.loader = Loader(chunk_size=self.chunk_size,
                             dry_run=self.dry_run)

    def close(self):
        super(ScoreProducer, self).close()
        self.loader.close()

    def process(self, data):
        target, disease, evidence, is_direct = data
        if evidence:
            score = self.scorer.score(target, disease, evidence, is_direct)
            if score.get_scoring_method(
                    ScoringMethods.HARMONIC_SUM).overall != 0:  # skip associations only with data with score 0

                # look for the gene in the lru cache
                gene_data = None
                hpa_data = None
                if target in self.lru_cache:
                    gene_data, hpa_data = self.lru_cache[target]

                if not gene_data:
                    gene_data = Gene()
                    try:
                        gene_data.load_json(
                            self.lookup_data.available_genes.get_gene(target,
                                                                      self.r_server))

                    except KeyError as e:
                        self.logger.debug('Cannot find gene code "%s" '
                                          'in lookup table' % target)
                        self.logger.exception(e)

                score.set_target_data(gene_data)

                # create a hpa expression empty jsonserializable class
                # to fill from Redis cache lookup_data
                if not hpa_data:
                    hpa_data = HPAExpression()
                    try:
                        hpa_data.update(
                            self.lookup_data.available_hpa.get_hpa(target,
                                                                   self.r_server))
                    except KeyError:
                        pass
                    except Exception as e:
                        self.logger.exception(e)

                    # set everything in the lru_cache
                    self.lru_cache[target] = (gene_data, hpa_data)

                try:
                    score.set_hpa_data(hpa_data)
                except KeyError:
                    pass
                except Exception as e:
                    self.logger.exception(e)

                disease_data = EFO()
                try:
                    disease_data.load_json(
                        self.lookup_data.available_efos.get_efo(disease, self.r_server))

                except KeyError as e:
                    self.logger.debug('Cannot find EFO code "%s" '
                                      'in lookup table' % disease)
                    self.logger.exception(e)

                score.set_disease_data(disease_data)

                if score: #bypass associations with overall score=0
                    element_id = '%s-%s' % (target, disease)
                    self.loader.put(Config.ELASTICSEARCH_DATA_ASSOCIATION_INDEX_NAME,
                                           Config.ELASTICSEARCH_DATA_ASSOCIATION_DOC_NAME,
                                           element_id,
                                           score,
                                           create_index=False)
                                           # routing=score.target['id']

                else:
                    logger.warning('Skipped association with score 0: %s-%s' % (target, disease))



class ScoreStorerWorker(RedisQueueWorkerProcess):
    def __init__(self,
                 score_q,
                 r_path,
                 chunk_size = 1e4,
                 dry_run = False,
                 es = None
                 ):
        super(ScoreStorerWorker, self).__init__(score_q, r_path)
        self.q = score_q
        self.chunk_size = chunk_size
        self.dry_run = dry_run
        self.es = None
        self.loader = None


    def process(self, data):

        target, disease, score = data
        element_id = '%s-%s' % (target, disease)
        self.loader.put(Config.ELASTICSEARCH_DATA_ASSOCIATION_INDEX_NAME,
                               Config.ELASTICSEARCH_DATA_ASSOCIATION_DOC_NAME,
                               element_id,
                               score,
                               create_index=False,
                               # routing=score.target['id'],
                            )

    def init(self):
        super(ScoreStorerWorker, self).init()
        self.loader = Loader(chunk_size=self.chunk_size,
                             dry_run=self.dry_run)

    def close(self):
        super(ScoreStorerWorker, self).close()
        self.loader.close()


class ScoringProcess():

    def __init__(self,
                 loader,
                 r_server):
        self.es_loader = loader
        self.es = loader.es
        self.es_query = ESQuery(loader.es)
        self.r_server = r_server

    def process_all(self,
                    targets = [],
                    dry_run = False):
        self.score_target_disease_pairs(targets=targets,
                                        dry_run=dry_run)

    def score_target_disease_pairs(self,
                                   targets = [],
                                   dry_run = False):

        overwrite_indices = not dry_run
        if overwrite_indices:
            overwrite_indices = not bool(targets)

        if not targets:
            targets =  list(self.es_query.get_all_target_ids_with_evidence_data())

        lookup_data = LookUpDataRetriever(self.es,
                                          self.r_server,
                                          targets=targets,
                                          data_types=(
                                              LookUpDataType.DISEASE,
                                              LookUpDataType.TARGET,
                                              LookUpDataType.ECO,
                                              LookUpDataType.HPA
                                          ),
                                          autoload=True if not targets
                                          or len(targets) > 100
                                          else False).lookup



        '''create queues'''
        number_of_workers = Config.WORKERS_NUMBER
        # too many storers
        number_of_storers = min(16, number_of_workers)
        queue_per_worker = 250
        if targets and len(targets) < number_of_workers:
            number_of_workers = len(targets)
        target_q = RedisQueue(queue_id=Config.UNIQUE_RUN_ID + '|target_q',
                              max_size=len(targets),
                              job_timeout=60*60*24,
                              r_server=self.r_server,
                              serialiser=None,
                              total=len(targets))
        disease_q = RedisQueue(queue_id=Config.UNIQUE_RUN_ID + '|disease_q',
                              max_size=len(lookup_data.diseases_with_data),
                              job_timeout=60*60*24,
                              r_server=self.r_server,
                              serialiser=None,
                              )

        target_map_q = RedisQueue(queue_id=Config.UNIQUE_RUN_ID + '|target_map_q',
                              max_size=len(targets),
                              job_timeout=60*60*24,
                              r_server=self.r_server,
                              serialiser='jsonpickle',
                              total=len(targets))

        target_disease_pair_q = RedisQueue(queue_id=Config.UNIQUE_RUN_ID + '|target_disease_pair_q',
                                           max_size=queue_per_worker * number_of_storers,
                                           job_timeout=1200,
                                           batch_size=10,
                                           r_server=self.r_server,
                                           serialiser='jsonpickle')
        # score_data_q = RedisQueue(queue_id=Config.UNIQUE_RUN_ID + '|score_data_q',
        #                           max_size=queue_per_worker * number_of_storers,
        #                           job_timeout=1200,
        #                           batch_size=10,
        #                           r_server=self.r_server,
        #                           serialiser='jsonpickle')

        q_reporter = RedisQueueStatusReporter([target_q,
                                               target_disease_pair_q,
                                               target_map_q,
                                               disease_q,
                                               # score_data_q
                                               ],
                                              interval=30,
                                              )
        q_reporter.start()


#         '''create data storage workers'''
#         storers = [ScoreStorerWorker(score_data_q,
#                                      None,
#                                      chunk_size=1000,
#                                      dry_run = dry_run,
#                                      ) for _ in range(number_of_storers)]
#
#         for w in storers:
#             w.start()

        # storage is located inside this code because the serialisation time
        scorers = [ScoreProducer(target_disease_pair_q,
                                 None,
                                 None,
                                 lookup_data,
                                 chunk_size=1000,
                                 dry_run=dry_run
                                 ) for _ in range(number_of_storers)]
        for w in scorers:
            w.start()


        '''start target-disease evidence producer'''
        readers = [TargetDiseaseEvidenceProducer(target_q,
                                                 None,
                                                 target_disease_pair_q,
                                                ) for _ in range(number_of_workers)]
        for w in readers:
            w.start()

        '''start  maps storers'''
        tm_storer = TargetMapStorer(target_map_q,
                                    None,
                                    disease_q,
                                    lookup_data=lookup_data,
                                    dry_run=dry_run
                                    )
        tm_storer.start()
        dm_storer = DiseaseMapStorer(disease_q,
                                    None,
                                    None,
                                    lookup_data=lookup_data,
                                    dry_run=dry_run
                                    )
        dm_storer.start()


        self.es_loader.create_new_index(Config.ELASTICSEARCH_DATA_ASSOCIATION_INDEX_NAME, recreate=overwrite_indices)
        self.es_loader.prepare_for_bulk_indexing(self.es_loader.get_versioned_index(Config.ELASTICSEARCH_DATA_ASSOCIATION_INDEX_NAME))



        for target in tqdm(targets,
                           desc='fetching evidence for targets',
                           unit=' targets',
                           file=tqdm_out,
                           unit_scale=True):
            target_q.put(target)
        target_q.set_submission_finished()

        '''wait for all workers to finish'''
        for w in readers:
            w.join()
        for w in scorers:
            w.join()
        tm_storer.join()
        dm_storer.join()

#         for w in storers:
#             w.join()

        logger.info('flushing data to index')
        self.es_loader.es.indices.flush('%s*'%Loader.get_versioned_index(Config.ELASTICSEARCH_DATA_ASSOCIATION_INDEX_NAME),
                                        wait_if_ongoing =True)


        logger.info("DONE")
