from yapsy.IPlugin import IPlugin
from mrtarget.modules.GeneData import Gene
from mrtarget.Settings import Config
from mrtarget.common.ElasticsearchQuery import ESQuery
from elasticsearch.exceptions import NotFoundError
from tqdm import tqdm
import sys
from mrtarget.common.LookupTables import _iterate_lut_file
import logging
logging.basicConfig(level=logging.INFO)

class Ensembl(IPlugin):

    def __init__(self, *args, **kwargs):
        self._logger = logging.getLogger(__name__)

    def print_name(self):
        self._logger.info("ENSEMBL gene data plugin")

    def merge_data(self, genes, loader, r_server, tqdm_out):
        for row in _iterate_lut_file(Config.ELASTICSEARCH_ENSEMBL_DOC_NAME):
            if row['id'] in genes:
                gene = genes.get_gene(row['id'])
                gene.load_ensembl_data(row)
                genes.add_gene(gene)
            else:
                gene = Gene()
                gene.load_ensembl_data(row)
                genes.add_gene(gene)

        self._clean_non_reference_genes(genes)

        self._logger.info("STATS AFTER ENSEMBL PARSING:\n" + genes.get_stats())

    def _clean_non_reference_genes(self, genes):
        for geneid, gene in genes.iterate():
            if not gene.is_ensembl_reference:
                genes.remove_gene(geneid)
