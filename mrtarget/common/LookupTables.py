import logging
from tqdm import tqdm
from mrtarget.common import TqdmToLogger
from mrtarget.common.ElasticsearchQuery import ESQuery
from mrtarget.common.Redis import RedisLookupTablePickle
from mrtarget.common.connection import new_redis_client, PipelineConnectors
from mrtarget.Settings import Config
from mrtarget.common import TqdmToLogger
import ujson as json
import glob
import os


def _build_lut_name(index_doc_name):
    '''return a list of data files with names'''
    return glob.glob(Config.OUTPUT_DIR + os.path.sep + \
                   index_doc_name + '_idx_data_' + Config.RELEASE_VERSION + '_*')


def _iterate_lut_file(index_doc_name):
    for filename in _build_lut_name(index_doc_name):
        with open(filename,'r') as f:
            for el in f:
                yield json.loads(el)


class HPALookUpTable(object):
    """
    A redis-based pickable hpa look up table using gene id as table
    id
    """

    def __init__(self,
                 es=None,
                 namespace=None,
                 r_server=None,
                 ttl=(60 * 60 * 24 + 7)):
        self.r_server = r_server
        self._table = RedisLookupTablePickle(namespace=namespace,
                                             r_server=self.r_server,
                                             ttl=ttl)
        self._logger = logging.getLogger(__name__)
        self.tqdm_out = TqdmToLogger(self._logger, level=logging.INFO)

        if self.r_server:
            self._load_hpa_data(self.r_server)

    def _load_hpa_data(self, r_server=None):
        self._logger.info("LUT expression loading data into redis...")

        for el in _iterate_lut_file(Config.ELASTICSEARCH_EXPRESSION_DOC_NAME):
            self.set_hpa(el, r_server=self._get_r_server(r_server))

        self._logger.info("LUT expression loading data into redis done")

    def get_hpa(self, idx, r_server=None):
        return self._table.get(idx, r_server=self._get_r_server(r_server))

    def set_hpa(self, hpa, r_server=None):
        self._table.set(hpa['gene'], hpa,
                        r_server=self._get_r_server(r_server))

    def get_available_hpa_ids(self, r_server=None):
        return self._table.keys(self._get_r_server(r_server))

    def __contains__(self, key, r_server=None):
        return self._table.__contains__(key,
                                        r_server=self._get_r_server(r_server))

    def __getitem__(self, key, r_server=None):
        return self.get_hpa(key, r_server=self._get_r_server(r_server))

    def __setitem__(self, key, value, r_server=None):
        self._table.set(key, value, r_server=self._get_r_server(r_server))

    def keys(self, r_server=None):
        return self._table.keys(self._get_r_server(r_server))

    def _get_r_server(self, r_server=None):
        return r_server if r_server else self.r_server


class EnsemblLookUpTable(object):
    """
    A redis-based pickable gene look up table
    """

    def __init__(self,
                 es=None,
                 namespace=None,
                 r_server=None,
                 ttl=60 * 60 * 24 + 7,
                 targets=[],
                 autoload=True):
        self._logger = logging.getLogger(__name__)
        self.r_server = r_server
        self._table = RedisLookupTablePickle(namespace=namespace,
                                             r_server=self.r_server,
                                             ttl=ttl)
        self._logger = logging.getLogger(__name__)
        if self.r_server and autoload:
            self.load_gene_data(self.r_server, targets)

    def load_gene_data(self, r_server=None, targets=[]):
        self._logger.info("LUT ensembl loading data into redis...")

        for target in _iterate_lut_file(Config.ELASTICSEARCH_ENSEMBL_DOC_NAME):
            self._table.set(target['id'], target, r_server=self._get_r_server(
                r_server))  # TODO can be improved by sending elements in batches

        self._logger.info("LUT ensembl loading data into redis done")

    def get_gene(self, target_id, r_server=None):
        try:
            return self._table.get(target_id, r_server=self._get_r_server(r_server))
        except KeyError:
            self._logger.exception('Cannot retrieve target from redis')
            raise KeyError()

    def set_gene(self, target, r_server=None):
        self._table.set(target['id'], target, r_server=self._get_r_server(r_server))

    def get_available_gene_ids(self, r_server=None):
        return self._table.keys(r_server=self._get_r_server(r_server))

    def __contains__(self, key, r_server=None):
        redis_contain = self._table.__contains__(key, r_server=self._get_r_server(r_server))
        if redis_contain:
            return True
        else:
            return False

    def __getitem__(self, key, r_server=None):
        return self.get_gene(key, self._get_r_server(r_server))

    def __setitem__(self, key, value, r_server=None):
        self._table.set(key, value, self._get_r_server(r_server))

    def __missing__(self, key):
        print key

    def keys(self, r_server=None):
        return self._table.keys(self._get_r_server(r_server))

    def _get_r_server(self, r_server=None):
        return r_server if r_server else self.r_server


class GeneLookUpTable(object):
    """
    A redis-based pickable gene look up table
    """

    def __init__(self,
                 es=None,
                 namespace = None,
                 r_server = None,
                 ttl = 60*60*24+7,
                 targets = [],
                 autoload=True):
        self._logger = logging.getLogger(__name__)
        self.r_server = r_server
        self._table = RedisLookupTablePickle(namespace = namespace,
                                            r_server = self.r_server,
                                            ttl = ttl)
        self._logger = logging.getLogger(__name__)
        self.uniprot2ensembl = {}
        if self.r_server and autoload:
            self.load_gene_data(self.r_server, targets)

    def load_gene_data(self, r_server = None, targets = []):
        self._logger.info("LUT gene loading data into redis...")

        for target in _iterate_lut_file(Config.ELASTICSEARCH_GENE_NAME_DOC_NAME):
            self._table.set(target['id'],target, r_server=self._get_r_server(r_server))#TODO can be improved by sending elements in batches
            if target['uniprot_id']:
                self.uniprot2ensembl[target['uniprot_id']] = target['id']
            for accession in target['uniprot_accessions']:
                self.uniprot2ensembl[accession] = target['id']

        self._logger.info("LUT gene loading data into redis done")

    def get_gene(self, target_id, r_server = None):
        try:
            return self._table.get(target_id, r_server=self._get_r_server(r_server))
        except KeyError:
            self._logger.exception('Cannot retrieve target from redis')
            raise KeyError()

    def set_gene(self, target, r_server = None):
        self._table.set(target['id'],target, r_server=self._get_r_server(r_server))

    def get_available_gene_ids(self, r_server = None):
        return self._table.keys(r_server = self._get_r_server(r_server))

    def __contains__(self, key, r_server=None):
        redis_contain = self._table.__contains__(key, r_server=self._get_r_server(r_server))
        if redis_contain:
            return True
        else:
            return False

    def __getitem__(self, key, r_server = None):
        return self.get_gene(key, self._get_r_server(r_server))

    def __setitem__(self, key, value, r_server=None):
        self._table.set(key, value, self._get_r_server(r_server))

    def __missing__(self, key):
        print key

    def keys(self, r_server=None):
        return self._table.keys(self._get_r_server(r_server))

    def _get_r_server(self, r_server = None):
        return r_server if r_server else self.r_server


class ECOLookUpTable(object):
    """
    A redis-based pickable gene look up table
    """


    def __init__(self,
                 es,
                 namespace=None,
                 r_server=None,
                 ttl=60 * 60 * 24 + 7):
        self._table = RedisLookupTablePickle(namespace=namespace,
                                             r_server=r_server,
                                             ttl=ttl)
        self.r_server = r_server
        self._logger = logging.getLogger(__name__)

        if r_server is not None:
            self._load_eco_data(r_server)

    @staticmethod
    def get_ontology_code_from_url(url):
        return url.split('/')[-1]

    def _load_eco_data(self, r_server=None):
        self._logger.info("LUT eco loading data into redis...")

        for eco in _iterate_lut_file(Config.ELASTICSEARCH_ECO_DOC_NAME):
            self._table.set(self.get_ontology_code_from_url(eco['code']), eco,
                            r_server=self._get_r_server(r_server))  # TODO can be improved by sending elements in batches

        self._logger.info("LUT eco loading data into redis done")

    def get_eco(self, efo_id, r_server=None):
        return self._table.get(efo_id, r_server=self._get_r_server(r_server))


    def set_eco(self, eco, r_server=None):
        self._table.set(self.get_ontology_code_from_url(eco['code']), eco, r_server=self._get_r_server(r_server))


    def get_available_eco_ids(self, r_server=None):
        return self._table.keys(r_server=self._get_r_server(r_server))


    def __contains__(self, key, r_server=None):
        return self._table.__contains__(key, r_server=self._get_r_server(r_server))


    def __getitem__(self, key, r_server=None):
        return self.get_eco(key, r_server=self._get_r_server(r_server))


    def __setitem__(self, key, value, r_server=None):
        self._table.set(key, value, r_server=self._get_r_server(r_server))


    def _get_r_server(self, r_server=None):
        return r_server if r_server else self.r_server


    def keys(self, r_server=None):
        return self._table.keys(r_server=self._get_r_server(r_server))

class EFOLookUpTable(object):
    """
    A redis-based pickable efo look up table.
    Allows to grab the EFO saved in ES and load it up in memory/redis so that it can be accessed quickly from multiple processes, reducing memory usage by sharing.
    """

    def __init__(self,
                 es=None,
                 namespace=None,
                 r_server=None,
                 ttl = 60*60*24+7):
        self.r_server = r_server
        self._table = RedisLookupTablePickle(namespace = namespace,
                                            r_server = self.r_server,
                                            ttl = ttl)
        self._logger = logging.getLogger(__name__)

        if self.r_server is not None:
            self._load_efo_data(r_server)

    def _load_efo_data(self, r_server = None):
        self._logger.info("LUT efo loading data into redis...")

        for efo in _iterate_lut_file(Config.ELASTICSEARCH_EFO_LABEL_DOC_NAME):
            self.set_efo(efo, r_server=self._get_r_server(r_server))

        self._logger.info("LUT efo loading data into redis done")

    def get_efo(self, efo_id, r_server=None):
        return self._table.get(efo_id, r_server=self._get_r_server(r_server))

    def set_efo(self, efo, r_server=None):
        efo_key = efo['path_codes'][0][-1]
        self._table.set(efo_key,efo, r_server=self._get_r_server(r_server))

    def get_available_gefo_ids(self, r_server=None):
        return self._table.keys(r_server=self._get_r_server(r_server))

    def __contains__(self, key, r_server=None):
        return self._table.__contains__(key, r_server=self._get_r_server(r_server))

    def __getitem__(self, key, r_server=None):
        return self.get_efo(key, r_server=self._get_r_server(r_server))

    def __setitem__(self, key, value, r_server=None):
        self._table.set(key, value, r_server=self._get_r_server(r_server))

    def keys(self, r_server=None):
        return self._table.keys(r_server=self._get_r_server(r_server))

    def _get_r_server(self, r_server = None):
        return r_server if r_server else self.r_server

class MPLookUpTable(object):
    """
    A redis-based pickable mp look up table.
    Allows to grab the MP saved in ES and load it up in memory/redis so that it can be accessed quickly from multiple processes, reducing memory usage by sharing.
    """

    def __init__(self,
                 es=None,
                 namespace=None,
                 r_server=None,
                 ttl = 60*60*24+7,
                 autoload=True):
        self.r_server = r_server
        self._table = RedisLookupTablePickle(namespace = namespace,
                                            r_server = self.r_server,
                                            ttl = ttl)

        self._logger = logging.getLogger(__name__)
        if self.r_server is not None and autoload:
            self._load_mp_data(r_server)

    def _load_mp_data(self, r_server = None):
        self._logger.info("LUT mp loading data into redis...")

        for mp in _iterate_lut_file(Config.ELASTICSEARCH_MP_LABEL_DOC_NAME):
            self.set_mp(mp, r_server=self._get_r_server(r_server))

        self._logger.info("LUT mp loading data into redis done")

    def get_mp(self, mp_id, r_server=None):
        return self._table.get(mp_id, r_server=self._get_r_server(r_server))

    def set_mp(self, mp, r_server=None):
        mp_key = mp['path_codes'][0][-1]
        self._table.set(mp_key, mp, r_server=self._get_r_server(r_server))

    def get_available_mp_ids(self, r_server=None):
        return self._table.keys(r_server=self._get_r_server(r_server))

    def __contains__(self, key, r_server=None):
        return self._table.__contains__(key, r_server=self._get_r_server(r_server))

    def __getitem__(self, key, r_server=None):
        return self.get_mp(key, r_server=self._get_r_server(r_server))

    def __setitem__(self, key, value, r_server=None):
        self._table.set(key, value, r_server=self._get_r_server(r_server))

    def keys(self, r_server=None):
        return self._table.keys(r_server=self._get_r_server(r_server))

    def _get_r_server(self, r_server = None):
        return r_server if r_server else self.r_server

class HPOLookUpTable(object):
    """
    A redis-based pickable hpo look up table.
    Allows to grab the HPO saved in ES and load it up in memory/redis so that it can be accessed quickly from multiple processes, reducing memory usage by sharing.
    """

    def __init__(self,
                 es=None,
                 namespace=None,
                 r_server=None,
                 ttl = 60*60*24+7):
        self.r_server = r_server
        self._table = RedisLookupTablePickle(namespace = namespace,
                                            r_server = self.r_server,
                                            ttl = ttl)
        self._logger = logging.getLogger(__name__)
        if self.r_server is not None:
            self._load_hpo_data(r_server)

    def _load_hpo_data(self, r_server = None):
        self._logger.info("LUT hpo loading data into redis...")

        for hpo in _iterate_lut_file(Config.ELASTICSEARCH_HPO_LABEL_DOC_NAME):
            self.set_hpo(hpo, r_server=self._get_r_server(r_server))

        self._logger.info("LUT hpo loading data into redis done")

    def get_hpo(self, hpo_id, r_server=None):
        return self._table.get(hpo_id, r_server=self._get_r_server(r_server))

    def set_hpo(self, hpo, r_server=None):
        hpo_key = hpo['path_codes'][0][-1]
        self._table.set(hpo_key, hpo, r_server=self._get_r_server(r_server))

    def get_available_hpo_ids(self, r_server=None):
        return self._table.keys(r_server=self._get_r_server(r_server))

    def __contains__(self, key, r_server=None):
        return self._table.__contains__(key, r_server=self._get_r_server(r_server))

    def __getitem__(self, key, r_server=None):
        return self.get_hpo(key, r_server=self._get_r_server(r_server))

    def __setitem__(self, key, value, r_server=None):
        self._table.set(key, value, r_server=self._get_r_server(r_server))

    def keys(self, r_server=None):
        return self._table.keys(r_server=self._get_r_server(r_server))

    def _get_r_server(self, r_server = None):
        return r_server if r_server else self.r_server

