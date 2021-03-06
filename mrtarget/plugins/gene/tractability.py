from yapsy.IPlugin import IPlugin
from mrtarget.Settings import Config
from tqdm import tqdm
from itertools import compress

import traceback
import logging
logging.basicConfig(level=logging.DEBUG)



class CancerBiomarkers(IPlugin):

    # Initiate Tractability object
    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self.loader = None
        self.r_server = None
        self.esquery = None
        self.ensembl_current = {}
        self.symbols = {}
        self.tractability = {}
        self.tqdm_out = None

    def print_name(self):
        self._logger.info("Tractability plugin")

    def merge_data(self, genes, loader, r_server, tqdm_out):

        self.loader = loader
        self.r_server = r_server
        self.tqdm_out = tqdm_out

        try:
            # Parse cancer biomarker data into self.cancerbiomarkers
            self.build_json(filename=Config.TRACTABILITY_FILENAME)

            # Iterate through all genes and add cancer biomarkers data if gene symbol is present
            self._logger.info("Tractability data injection")
            for gene_id, gene in tqdm(genes.iterate(),
                                      desc='Adding Tractability data',
                                      unit=' gene',
                                      file=self.tqdm_out):
                # Extend gene with related Cancer Biomarker data
                if gene.ensembl_gene_id in self.tractability:
                    self._logger.debug("Adding Cancer Biomarker data to gene %s", gene.ensembl_gene_id)
                    gene.tractability=self.tractability[gene.ensembl_gene_id]

        except Exception as ex:
            self._logger.exception(str(ex), exc_info=1)
            raise ex

    def build_json(self, filename=Config.TRACTABILITY_FILENAME):

        sm_bucketList = [1, 2, 3, 4, 5, 6, 7, 8]
        ab_bucketList = [1, 2, 3, 4, 5, 6, 7, 8, 9]

        with open(filename, 'r') as input:
            next(input)
            for row in input:
                (ensembl_gene_id, accession, Bucket_1, Bucket_2, Bucket_3, Bucket_4, Bucket_5, Bucket_6, Bucket_7,
                 Bucket_8,
                 Bucket_sum, Top_bucket, Category, Clinical_Precedence, Discovery_Precedence, Predicted_Tractable,
                 ensemble, High_Quality_ChEMBL_compounds, Small_Molecule_Druggable_Genome_Member,
                 Bucket_1_ab, Bucket_2_ab, Bucket_3_ab, Bucket_4_ab, Bucket_5_ab, Bucket_6_ab, Bucket_7_ab, Bucket_8_ab,
                 Bucket_9_ab, Bucket_sum_ab, Top_bucket_ab, Uniprot_high_conf_loc, GO_high_conf_loc,
                 Uniprot_med_conf_loc,
                 GO_med_conf_loc, Transmembrane, Signal_peptide, HPA_main_location, Clinical_Precedence_ab,
                 Predicted_Tractable__High_confidence, Predicted_Tractable__Medium_to_low_confidence, Category_ab) = \
                    tuple(row.rstrip().split('\t'))

                # Get lists of small molecule and antibody buckets
                sm_buckets = list(compress(sm_bucketList, [x == '1' for x in
                                                        [Bucket_1, Bucket_2, Bucket_3, Bucket_4, Bucket_5, Bucket_6,
                                                         Bucket_7, Bucket_8]]))
                ab_buckets = list(compress(ab_bucketList, [x == '1' for x in
                                                        [Bucket_1_ab, Bucket_2_ab, Bucket_3_ab, Bucket_4_ab,
                                                         Bucket_5_ab, Bucket_6_ab, Bucket_7_ab, Bucket_8_ab,
                                                         Bucket_9_ab]]))

                line = {'smallmolecule': {}, 'antibody': {}}
                line['smallmolecule'] = {
                    'buckets': sm_buckets,  # list of buckets
                    'categories': {
                        'clinical_precedence': Clinical_Precedence,
                        'discovery_precedence': Discovery_Precedence,
                        'predicted_precedence': Predicted_Tractable
                    },
                    'top_category': Category,
                    'ensemble': ensemble,
                    'high_quality_compounds': High_Quality_ChEMBL_compounds,
                    'small_molecule_genome_member': Small_Molecule_Druggable_Genome_Member
                }
                line['antibody'] = {
                    'buckets': ab_buckets,
                    'categories': {
                        'clinical_precedence': Clinical_Precedence_ab,
                        'predicted_tractable_high_confidence': Predicted_Tractable__High_confidence,
                        'predicted_tractable_med_low_confidence': Predicted_Tractable__Medium_to_low_confidence
                    },
                    'top_category': Category_ab
                }

                # Add data for current biomarker to self.cancerbiomarkers
                self.tractability[ensembl_gene_id] = line