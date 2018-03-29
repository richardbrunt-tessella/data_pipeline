import copy
import re
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import pysftp
import gzip
import pickle
from paramiko import AuthenticationException
import logging
import ujson
import rdflib
import requests
from rdflib import URIRef
from rdflib.namespace import Namespace, NamespaceManager
from rdflib.namespace import OWL, RDF, RDFS
from mrtarget.common import Actions
from SPARQLWrapper import SPARQLWrapper, JSON
from tqdm import tqdm
from mrtarget.common import TqdmToLogger
from datetime import date
from mrtarget.Settings import Config
from ontologyutils import rdf_utils


__copyright__ = "Copyright 2014-2018, Open Targets"
__credits__ = []
__license__ = "Apache 2.0"
__version__ = ""
__maintainer__ = "Gautier Koscielny"
__email__ = "gautier.x.koscielny@gsk.com"
__status__ = "Production"


logger = logging.getLogger(__name__)


class OntologyActions(Actions):
    CHECKDISEASES = 'checkdiseases'
    PHENOTYPESLIM = 'phenotypeslim'
    DISEASEPHENOTYPES = 'diseasephenotypes'



class DiseasePhenotypeReader():

    def __init__(self):
        self.lookup_data = None
        self.efo_ontology = None
        self.hpo_ontology = None
        self.mp_ontology = None
        self.invalid_diseases = dict()
        self.obsolete_diseases = dict()

        self.logger = logging.getLogger(__name__)

    def get_ontologies(self):

        self.efo_ontology = rdf_utils.OntologyClassReader().get_efo()
        self.hpo_ontology = rdf_utils.OntologyClassReader().get_hpo()
        self.mp_ontology = rdf_utils.OntologyClassReader().get_mp()

        # print the first 10 classes from EFO
        c = 0
        for k,v in self.efo_ontology.current_classes.iteritems():
            print(k)
            c+=1
            if c > 10:
                break

        if 'http://www.ebi.ac.uk/efo/EFO_0003869' in self.efo_ontology.current_classes:
            print('found breast neoplasm')
        else:
            print('NO http://www.ebi.ac.uk/efo/EFO_0003869')

    def parse_gzipfile(self, file_path):

        file_size, file_mod_time = os.path.getsize(file_path), os.path.getmtime(file_path)
        with open(file_path, mode='rb') as f:
            with gzip.GzipFile(filename=file_path,
                               mode='rb',
                               fileobj=f,
                               mtime=file_mod_time) as fh:
    
                self.logger.info('Starting parsing %s' % file_path)
    
                line_buffer = []
                offset = 0
                chunk = 1
                line_number = 0
    
                for line in fh:
                    python_raw = ujson.loads(line)
                    # get disease IRI
                    disease_iri = python_raw['disease']['id']
                    # check that the disease id is in EFO
                    self.check_disease_in_evidence(disease_iri, file_path)
                    #if re.match('http://purl.obolibrary.org/obo/HP_\d+', id) or re.match(
                    #        'http://purl.obolibrary.org/obo/MP_\d+', id):
                    #    ''' get all terms '''
                    #    self.phenotype_map[id] = self.phenotypes.classes_paths[id]

        for k,v in self.invalid_diseases.items():
            print("invalid disease %s from %s"%(k, '; '.join(list(v['sources']))))
        for k,v in self.obsolete_diseases.items():
            print("obsolete disease %s from %s"%(k, '; '.join(list(v['sources']))))

    def check_disease_in_evidence(self, disease_iri, file_path):

        # Check disease term or phenotype term
        if (disease_iri not in self.efo_ontology.current_classes) and \
                (disease_iri not in self.hpo_ontology.current_classes) and \
                (disease_iri not in self.mp_ontology.current_classes):  # or \
            # (disease_id in self.efo_uncat):
            if disease_iri not in self.invalid_diseases:
                self.invalid_diseases[disease_iri] = dict(sources=set(file_path))
            else:
                self.invalid_diseases[disease_iri]['sources'].add(file_path)

        if (disease_iri in self.efo_ontology.obsolete_classes) or \
                (disease_iri in self.hpo_ontology.obsolete_classes) or \
                (disease_iri in self.mp_ontology.obsolete_classes):
            if disease_iri not in self.obsolete_diseases:
                self.obsolete_diseases[disease_iri] = dict(sources=set(file_path))
            else:
                self.obsolete_diseases[disease_iri]['sources'].add(file_path)

class PhenotypeSlim():

    def __init__(self):

        self.rdf_graph = rdflib.Graph()
        self.phenotypes = None
        #self.efo = OntologyClassReader()
        #self.efo.load_efo_classes()

        self.phenotype_map = {}
        self.phenotype_excluded = set()

        self._remote_filenames = dict()
        self._logger = logging.getLogger(self.__class__.__name__)
        tqdm_out = TqdmToLogger(self._logger, level=logging.INFO)


    def get_ontology_path(self, base_class, term):

        """ if the term is already there or if it is a disease term """
        if term in self.phenotype_map:
            return


        #if term == 'http://purl.obolibrary.org/obo/HP_0001251':
        if True:

            for sparql_query in [DIRECT_ANCESTORS, INDIRECT_ANCESTORS]:
                self.sparql.setQuery(sparql_query%(term, term))
                self.sparql.setReturnFormat(JSON)
                results = None
                n = 0
                while (n < 2):
                    try:
                        results = self.sparql.query().convert()
                        n = 3
                    except SPARQLWrapper.SPARQLExceptions.EndPointNotFound, e:
                        print e
                        self._logger.error(e)
                        if n > 2:
                            raise e
                        else:
                            n=n+1

                #print len(results)
                #print json.dumps(results)

                for result in results["results"]["bindings"]:
                    #print json.dumps(result)
                    count = int(result['distance']['value'])
                    parent_label = result['ancestor_label']['value']
                    ancestor = result['ancestor']['value']
                    direct_child = result['direct_child']['value']
                    direct_child_label = result['direct_child_label']['value']
                    ''' Add the term to the phenotype map '''
                    if direct_child not in self.phenotype_map:
                        self.phenotype_map[direct_child] = { 'label': direct_child_label , 'superclasses': [] }
                    ''' Put all the ancestors to the phenotype map '''
                    if ancestor not in self.phenotype_map[direct_child]['superclasses']:
                        self.phenotype_map[direct_child]['superclasses'].append(ancestor)
                        print "%i %s %s (direct child is %s %s)"%(count, parent_label, ancestor, direct_child_label, direct_child)
                        print "---------"
                    #print "%i %s %s (direct child is %s %s)"%(count, parent_label, ancestor, direct_child_label, direct_child)


    def load_all_phenotypes(self):
        '''
        Load HPO and MP to accept phenotype terms that are not in EFO
        :return:
        '''
        self.phenotypes = OntologyClassReader()
        self.phenotypes.load_hpo_classes()
        self.phenotypes.get_classes_paths(root_uri='http://purl.obolibrary.org/obo/HP_0000118')
        self.phenotypes.load_mp_classes()
        self.phenotypes.get_classes_paths(root_uri='http://purl.obolibrary.org/obo/MP_0000001')

    def exclude_phenotypes(self, l):
        '''
        :param l:
        :return:
        '''
        for p in l:
            if p not in self.phenotype_excluded:
                self.phenotype_excluded.add(p)
                print "Excluding %s"%p
                # get parents
                sparql_query = DIRECT_ANCESTORS
                self.sparql.setQuery(sparql_query%(p, p))
                self.sparql.setReturnFormat(JSON)
                results = self.sparql.query().convert()
                al = []
                for result in results["results"]["bindings"]:
                    count = int(result['distance']['value'])
                    parent_label = result['ancestor_label']['value']
                    ancestor = result['ancestor']['value']
                    al.append(ancestor)
                    self.exclude_phenotypes(al)

    def _store_remote_filename(self, filename):
        # print "%s" % filename
        self._logger.debug("%s" % filename)
        if filename.startswith('/upload/submissions/') and \
            filename.endswith('.json.gz'):
            self._logger.debug("%s" % filename)
            if True:
                version_name = filename.split('/')[3].split('.')[0]
                # print "%s" % filename
                if '-' in version_name:
                    user, day, month, year = version_name.split('-')
                    if '_' in user:
                        datasource = ''.join(user.split('_')[1:])
                        user = user.split('_')[0]
                    else:
                        datasource = Config.DATASOURCE_INTERNAL_NAME_TRANSLATION_REVERSED[user]
                    release_date = date(int(year), int(month), int(day))

                    if user not in self._remote_filenames:
                        self._remote_filenames[user]={datasource : dict(date = release_date,
                                                                          file_path = filename,
                                                                          file_version = version_name)
                                                      }
                    elif datasource not in self._remote_filenames[user]:
                        self._remote_filenames[user][datasource] = dict(date=release_date,
                                                                        file_path=filename,
                                                                        file_version=version_name)
                    else:
                        if release_date > self._remote_filenames[user][datasource]['date']:
                            self._remote_filenames[user][datasource] = dict(date=release_date,
                                                                file_path=filename,
                                                                file_version=version_name)
            #except e:
            #    self.logger.error("%s Error checking file %s: %s" % (self.__class__.__name__, filename, e))
            #    print 'error getting remote file%s'%filename
            #    self.logger.debug('error getting remote file%s'%filename)

    def _callback_not_used(self, path):
        self._logger.debug("skipped " + path)

    def create_phenotype_slim_from_selection(self):

        for url in Config.PHENOTYPE_SLIM_INPUT_URLS:
            response = requests.get(url)
            self._logger.info("Read url %s - response code %s" % (url, response.code))
            lines = response.readlines()

            for line in lines:
                self._logger.info(line.rstrip())

    def create_phenotype_slim_from_evidence(self, local_files = []):

        self.load_all_phenotypes()

        if local_files:

            for file_path in local_files:
                self._logger.info("Parsing file %s" % (file_path))
                file_size, file_mod_time = os.path.getsize(file_path), os.path.getmtime(file_path)
                with open(file_path, mode='rb') as f:
                    self.parse_gzipfile(filename=file_path, mode='rb', fileobj=f, mtime=file_mod_time)
        else:

            for u in tqdm(Config.ONTOLOGY_PREPROCESSING_FTP_ACCOUNTS,
                             desc='scanning ftp accounts',
                             # file=tqdm_out,
                             leave=False):
                try:
                    p = Config.EVIDENCEVALIDATION_FTP_ACCOUNTS[u]
                    self._logger.info("%s %s" % (u, p))
                    cnopts = pysftp.CnOpts()
                    cnopts.hostkeys = None  # disable host key checking.
                    with pysftp.Connection(host=Config.EVIDENCEVALIDATION_FTP_HOST['host'],
                                           port=Config.EVIDENCEVALIDATION_FTP_HOST['port'],
                                           username=u,
                                           password=p,
                                           cnopts = cnopts,
                                           ) as srv:
                        srv.walktree('/', fcallback=self._store_remote_filename, dcallback=self._callback_not_used, ucallback=self._callback_not_used)
                        srv.close()
                        for datasource, file_data in tqdm(self._remote_filenames[u].items(),
                                                          desc='scanning available datasource for account %s'%u,
                                                          # file=tqdm_out,
                                                          leave=False,):
                            latest_file = file_data['file_path']
                            file_version = file_data['file_version']
                            self._logger.info("found latest file %s for datasource %s" % (latest_file, datasource))
                            self.parse_gzipfile(latest_file, u, p)
                except AuthenticationException:
                    self._logger.error('cannot connect with credentials: user:%s password:%s' % (u, p))

        for uri, p in self.phenotype_map.iteritems():
            logger.debug(uri)
            logger.debug(json.dumps(p, indent=2))




def main():

    obj = DiseasePhenotypeReader()
    obj.get_ontologies()
    obj.parse_gzipfile(file_path='/Users/otvisitor/Documents/data/18.02/phewas_catalog-11-09-2017.json.gz')
    #obj.parse_gzipfile(file_path='/Users/otvisitor/Documents/data/18.02/chembl-02-02-2018.json.gz')

    #obj = PhenotypeSlim()
    #obj.load_all_phenotypes()
    #obj.create_phenotype_slim(local_files=['/Users/koscieln/Documents/work/gitlab/data_pipeline/samples/cttv008-22-07-2016.json.gz'])


if __name__ == "__main__":
    main()
