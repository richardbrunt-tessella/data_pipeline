from __future__ import print_function
import logging
import argparse
import sys
import itertools as it
from logging.config import fileConfig

from mrtarget.common.Redis import enable_profiling
from mrtarget.common.ElasticsearchLoader import Loader
from mrtarget.common.connection import PipelineConnectors
from mrtarget.ElasticsearchConfig import ElasticSearchConfiguration
from mrtarget.modules.Association import ScoringProcess
from mrtarget.modules.DataDrivenRelation import DataDrivenRelationProcess
from mrtarget.modules.Dump import DumpGenerator
from mrtarget.modules.ECO import EcoProcess
from mrtarget.modules.EFO import EfoProcess
from mrtarget.modules.HPO import HpoProcess
from mrtarget.modules.MP import MpProcess
from mrtarget.modules.Ensembl import EnsemblProcess
from mrtarget.modules.EvidenceString import EvidenceStringProcess
from mrtarget.modules.EvidenceValidation import EvidenceValidationFileChecker
from mrtarget.modules.GeneData import GeneManager
from mrtarget.modules.HPA import HPAProcess
from mrtarget.modules.MouseModels import MouseModelsActions, Phenodigm
from mrtarget.modules.QC import QCRunner
from mrtarget.modules.Reactome import ReactomeProcess
from mrtarget.modules.SearchObjects import SearchObjectProcess
from mrtarget.modules.Uniprot import UniprotDownloader
from mrtarget.modules.Metrics import Metrics
from mrtarget.Settings import Config, file_or_resource, update_schema_version

def main():

    #set up logging
    fileConfig(file_or_resource('logging.ini'),  disable_existing_loggers=False)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='Open Targets processing pipeline')

    #take the release tag from the command line, but fall back to environment or ini files
    parser.add_argument('release_tag', nargs='?', default=Config.RELEASE_VERSION,
                        help='The prefix to prepend default: %s' % \
                        Config.RELEASE_VERSION)

    #load supplemental and genetic informtaion from various external resources
    parser.add_argument("--hpa", help="download human protein atlas, process, and store in elasticsearch",
                        action="store_true")
    parser.add_argument("--ens", help="retrieve the latest ensembl gene records, store in elasticsearch",
                        action="store_true")
    parser.add_argument("--unic", help="cache the uniprot human entries in elasticsearch",
                        action="store_true")
    parser.add_argument("--rea", help="download reactome data, process it, and store elasticsearch",
                        action="store_true")

    #use the sources to combine the gene information into a single new index
    parser.add_argument("--gen", help="merge the available gene information, store in elasticsearch",
                        action="store_true")

    #load various ontologies into various indexes
    parser.add_argument("--mp", help="process Mammalian Phenotype (MP), store the resulting json objects in elasticsearch",
                         action="store_true")
    parser.add_argument("--efo", help="process Experimental Factor Ontology (EFO), store in elasticsearch",
                        action="store_true")
    parser.add_argument("--eco", help="process Evidence and Conclusion Ontology (ECO), store in elasticsearch",
                        action="store_true")
    parser.add_argument("--hpo", help="process Human Phenotype Ontology (HPO), store in elasticsearch",
                         action="store_true")

    #this generates a elasticsearch index from a source json file
    parser.add_argument("--val", help="check json file, validate, and store in elasticsearch",
                        action="store_true")
    parser.add_argument("--valreset", help="reset audit table and previously parsed evidencestrings",
                        action="store_true")
    parser.add_argument("--input-file", help="pass the path to a gzipped file to use as input for the data validation step",
                        action='append', default=[])
    parser.add_argument("--schema-version", help="set the schema version aka 'branch' name. Default is 'master'",
                        action='store', default='master')

    #this is related to generating a combine evidence index from all the inidividual datasource indicies
    parser.add_argument("--evs", help="process and validate the available evidence strings, store in elasticsearch",
                        action="store_true")
    parser.add_argument("--datasource", help="just process data for this datasource. Does not work with all the steps!!",
                        action='append', default=[])

    #this has to be stored as "ass" instead of "as" because "as" is a reserved name when accessing it later e.g. `args.as`
    parser.add_argument("--as", help="compute association scores, store in elasticsearch",
                        action="store_true", dest="ass")                        
    parser.add_argument("--targets", help="just process data for this target. Does not work with all the steps!!",
                        action='append', default=[])
                        
    #these are related to generated in a search index
    parser.add_argument("--sea", help="compute search results, store in elasticsearch",
                        action="store_true")
    parser.add_argument("--skip-diseases", help="Skip adding diseases to the search index",
                        action='store_true', default=False)
    parser.add_argument("--skip-targets", help="Skip adding targets to the search index",
                        action='store_true', default=False)

    #additional information to add
    parser.add_argument("--ddr", help="compute data driven t2t and d2d relations, store in elasticsearch",
                        action="store_true")

    #generate some high-level summary metrics over the release
    parser.add_argument("--metric", help="generate metrics", 
                        action="store_true")

    #quality control steps
    parser.add_argument("--qc", help="Run quality control scripts",
                        action="store_true")

    #phenodigm related processes
    parser.add_argument("--musu", dest='mus', help="update mouse model data",
                        action="append_const", const = MouseModelsActions.UPDATE_CACHE)
    parser.add_argument("--musg", dest='mus', help="update mus musculus and home sapiens gene list",
                        action="append_const", const = MouseModelsActions.UPDATE_GENES)
    parser.add_argument("--muse", dest='mus', help="generate mouse model evidence",
                        action="append_const", const = MouseModelsActions.GENERATE_EVIDENCE)
    parser.add_argument("--mus", dest='mus', help="update mouse models data",
                        action="append_const", const = MouseModelsActions.ALL)
                       
    #use an external redis rather than spawning one ourselves
    parser.add_argument("--persist-redis", help="the temporary file wont be deleted if True default: False",
                        action='store_true', default=False)
    parser.add_argument("--redis-remote", help="connect to a remote redis",
                        action='store_true', default=False)
    parser.add_argument("--redis-host", help="redis host",
                        action='store', default='')
    parser.add_argument("--redis-port", help="redis port",
                        action='store', default='')

    #tweak how lookup tables are managed
    parser.add_argument("--lt-reuse", help="reuse the current lookuptable",
                        action='store_true', default=False)
    parser.add_argument("--lt-namespace", help="lookuptable namespace to reuse",
                        action='store', default='')

    #for debugging
    parser.add_argument("--dump", help="dump core data to local gzipped files",
                        action="store_true")
    parser.add_argument("--dry-run", help="do not store data in the backend, useful for dev work. Does not work with all the steps!!",
                        action='store_true', default=False)
    parser.add_argument("--profile", help="magically profiling process() per process",
                        action='store_true', default=False)
    parser.add_argument("--log-level", help="set the log level",
                        action='store', default='WARNING')
                        
    args = parser.parse_args()

    if not args.release_tag:
        logger.error('A [release-tag] has to be specified.')
        print('A [release-tag] has to be specified.', file=sys.stderr)
        return 1
    else:
        Config.RELEASE_VERSION = args.release_tag

    targets = args.targets

    if args.lt_namespace:
        Config.LT_NAMESPACE = args.lt_namespace

    if args.lt_reuse:
        Config.LT_REUSE = True

    if args.redis_remote:
        Config.REDISLITE_REMOTE = args.redis_remote

    if args.redis_host:
        Config.REDISLITE_DB_HOST = args.redis_host

    if args.redis_port:
        Config.REDISLITE_DB_PORT = args.redis_port

    enable_profiling(args.profile)

    logger.debug('redis remote %s and host %s port %s',
                 str(Config.REDISLITE_REMOTE),
                 Config.REDISLITE_DB_HOST,
                 Config.REDISLITE_DB_PORT)

    connectors = PipelineConnectors()

    if args.log_level:
        try:
            root_logger = logging.getLogger()
            root_logger.setLevel(logging.getLevelName(args.log_level))
            logger.setLevel(logging.getLevelName(args.log_level))
            logger.info('main log level set to: '+ str(args.log_level))
            root_logger.info('root log level set to: '+ str(args.log_level))
        except Exception, e:
            root_logger.exception(e)
            return 1

    connected = connectors.init_services_connections(redispersist=args.persist_redis)

    logger.debug('Attempting to establish connection to the backend... %s',
                 str(connected))

    logger.info('setting release version %s' % Config.RELEASE_VERSION)


    with Loader(connectors.es,
                chunk_size=ElasticSearchConfiguration.bulk_load_chunk,
                dry_run = args.dry_run) as loader:

        # get the schema version and change all needed resources
        update_schema_version(Config,args.schema_version)
        logger.info('setting schema version string to %s', args.schema_version)

        if args.rea:
            ReactomeProcess(loader).process_all()
        if args.ens:
            EnsemblProcess(loader).process()
        if args.unic:
            UniprotDownloader(loader).cache_human_entries()
        if args.hpa:
            HPAProcess(loader,connectors.r_server).process_all(dry_run=args.dry_run)

        if args.gen:
            GeneManager(loader,connectors.r_server).merge_all(dry_run=args.dry_run)

        if args.mp:
            MpProcess(loader).process_all()
        if args.efo:
            EfoProcess(loader).process_all()
        if args.eco:
            EcoProcess(loader).process_all()
        if args.hpo:
            HpoProcess(loader).process_all()

        if args.mus:
            do_all = (MouseModelsActions.ALL in args.mus)
            if (MouseModelsActions.UPDATE_CACHE in args.mus) or do_all:
                Phenodigm(connectors.es, connectors.r_server).update_cache()
            if (MouseModelsActions.UPDATE_GENES in args.mus) or do_all:
                Phenodigm(connectors.es, connectors.r_server).update_genes()
            if (MouseModelsActions.GENERATE_EVIDENCE in args.mus) or do_all:
                Phenodigm(connectors.es, connectors.r_server).generate_evidence()


        if args.val:
            if args.input_file:
                input_files = list(it.chain.from_iterable([el.split(",") for el in args.input_file]))
            else:
                #default behaviour: use all the data sources listed in the evidences_sources.txt file
                logger.debug('reading the evidences sources URLs from evidence_sources.txt')
                with open(file_or_resource('evidences_sources.txt')) as f:
                    input_files = [x.rstrip() for x in f.readlines()]
            EvidenceValidationFileChecker(connectors.es, connectors.r_server, 
                dry_run=args.dry_run).check_all(input_files=input_files)
        if args.valreset:
            EvidenceValidationFileChecker(connectors.es, connectors.r_server).reset()

        if args.evs:
            targets = EvidenceStringProcess(connectors.es,
                                                connectors.r_server,
                                                ).process_all(datasources = args.datasource,
                                                                                      dry_run=args.dry_run)
        if args.ass:
            ScoringProcess(loader, connectors.r_server).process_all(targets = targets,
                                                             dry_run=args.dry_run)
        if args.ddr:
            DataDrivenRelationProcess(connectors.es, connectors.r_server).process_all(dry_run = args.dry_run)

        if args.sea:
            SearchObjectProcess(loader, connectors.r_server).process_all(skip_targets=args.skip_targets,
                                                                             skip_diseases=args.skip_diseases)
        if args.metric:
            Metrics(connectors.es).generate_metrics()

        if args.qc:
            QCRunner(connectors.es).run_associationQC()

        if args.dump:
            DumpGenerator().dump()

    logger.debug('close connectors')
    connectors.close()

    logger.info('`'+" ".join(sys.argv)+'` - finished')
    return 0


if __name__ == '__main__':
    sys.exit(main())
