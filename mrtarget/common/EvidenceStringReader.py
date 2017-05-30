import logging
import json
from addict import Dict
from mrtarget.Settings import Config
from mrtarget.common import generate_validator_from_schema, URLZSource


__copyright__ = "Copyright 2014-2016, Open Targets"
__credits__ = []
__license__ = "Apache 2.0"
__version__ = ""
__maintainer__ = "Gautier Koscielny"
__email__ = "gautierk@opentargets.org"
__status__ = "Production"

from logging.config import fileConfig

logger = logging.getLogger('root')


class EvidenceStringReader(object):
    def __init__(self):
        pass

    def parse_gzipfile(self, filename, out_file):

        out_fh = open(out_file, 'w')

        with URLZSource(filename).open() as fh:

            logging.info('Starting parsing %s' % filename)

            for line in fh:
                out_line = line.rstrip()
                python_raw = json.loads(line)
                obj = None
                data_type = python_raw['type']
                if data_type in Config.EVIDENCEVALIDATION_DATATYPES:
                    if data_type == 'genetic_association':
                        uri = Config.EVIDENCEVALIDATION_VALIDATOR_SCHEMAS[data_type]
                        validator = generate_validator_from_schema(uri)
                        validation_errors = [str(e) for e in \
                                             validator.iter_errors(python_raw)]

                        # check evidence after is validation_errors
                        if not validation_errors:
                            obj = Dict(python_raw)
                            if obj.evidence.gene2variant.resource_score is None:
                                obj.evidence.gene2variant.resource_score.type = "probability"
                                obj.evidence.gene2variant.resource_score.method.description = "NA"
                                obj.evidence.gene2variant.resource_score.method.reference = "NA"
                                obj.evidence.gene2variant.resource_score.method.url = "NA"
                                obj.evidence.gene2variant.resource_score.vaue = 1.
                                # from addict to python dict
                                out_line = json.dumps(obj.to_dict())
                            else:
                                print obj.evidence.gene2variant.resource_score.value
                        else:
                            print "error: " + '\n'.join(validation_errors)
                out_fh.write(out_line + "\n")
        fh.close()
        out_fh.close()

def main():

    obj = EvidenceStringReader()
    obj.parse_gzipfile(filename='file:///Users/koscieln/Documents/data/ftp/cttv012/upload/submissions/cttv012-22-11-2016.json.gz', out_file='/Users/koscieln/Documents/data/ftp/cttv012/upload/submissions/cttv012-28-11-2016.json')


if __name__ == "__main__":
    main()
