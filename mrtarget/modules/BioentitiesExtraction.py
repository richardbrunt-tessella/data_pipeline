import sys
import httplib
import time
import optparse
from tqdm import tqdm
import requests
from mrtarget.common import TqdmToLogger
import logging
import os
import json
import re
import hashlib
import datetime



__copyright__ = "Copyright 2014-2017, Open Targets"
__credits__ = ["Gautier Koscielny"]
__license__ = "Apache 2.0"
__version__ = "1.2.6"
__maintainer__ = "Gautier Koscielny"
__email__ = "gautierk@opentargets.org"
__status__ = "Production"




class ProteinComplexExtractor(object):

    def __init__(self):
        self._logger = logging.getLogger(__name__)
        tqdm_out = TqdmToLogger(self._logger,level=logging.INFO)
        self.chembl = dict()
        self.corum = dict()
        self.go = dict()


    def parse_ChEMBL(self):

        # get number of entries
        total_count = 0
        query = '/chembl/api/data/target/?offset=0&limit=100'
        offset = 0
        while True:

            url = 'https://www.ebi.ac.uk%s&format=json'% (query)
            #self._logger.info(url)
            #print url
            r = requests.get(url)
            results = r.json()
            for target in results["targets"]:
                if target["organism"] == "Homo sapiens" and target["target_type"] in ["PROTEIN COMPLEX", "PROTEIN FAMILY", "PROTEIN COMPLEX GROUP"]:
                    target_chembl_id = target["target_chembl_id"]
                    pref_name = target["pref_name"]
                    self.chembl[target_chembl_id] = { 'name' : pref_name }
                    for target_component in target["target_components"]:
                        print "%s\t%s\t%s\t%s\t%s\t%s" %(target_chembl_id, pref_name, target["target_type"], target_component["component_type"], target_component["relationship"], target_component["accession"])
            #self._logger.info(results["page_meta"]["next"])

            if not results["page_meta"]["next"]:
                break
            else:
                query = results["page_meta"]["next"]

    def parse_CORUM(self):
        url = 'http://mips.helmholtz-muenchen.de/corum/download/allComplexes.json'
        r = requests.get(url)
        results = r.json()
        total = 0
        for complex in results:
            total+=1
            self.corum["%i"%(complex["ComplexID"])] = complex
        print total

    def parse_GO(self, next_url=None):

        if next_url:
            url = next_url
        else:
            term_iri = 'http://purl.obolibrary.org/obo/GO_0032991'
            print requests.utils.quote(term_iri, safe='')
            url = 'http://www.ebi.ac.uk/ols/api/ontologies/go/terms/%s'%(requests.utils.quote(requests.utils.quote(term_iri, safe='')))
        print url
        r = requests.get(url)
        results = r.json()
        terms = list()
        if "_embedded" in results:
            terms = results["_embedded"]["terms"]
        else:
            terms.append(results)
        for term in terms:
            obo_id = term["obo_id"].replace(':', '_')
            if obo_id not in self.go:
                print json.dumps(term, indent=2)
                obo_label = term["label"]
                self.go[obo_id] = [ obo_label ]
                if term["synonyms"]:
                    self.go[obo_id].extend(term["synonyms"])
                if "hierarchicalChildren" in term["_links"]:
                    self.parse_GO(next_url=term["_links"]["hierarchicalChildren"]["href"])



    def generate_bioentities_files(self):

        filename = "chembl"
        chembl_json = dict()
        for k,v in self.chembl.items():
            chembl_json[v["name"]] = [ k ]
        with open('/Users/otvisitor/Documents/data/ChEMBL_proteincomplexes.json', 'wb') as outfile:
            json.dump(chembl_json, outfile)

        corum_json = dict()
        for k,v in self.corum.items():
            corum_json[v["ComplexName"]] = [ k ]
            if v["Synonyms"]:
                synonyms = v["Synonyms"].split("; ")
                for synonym in synonyms:
                    corum_json[synonym] = [ k ]

        with open('/Users/otvisitor/Documents/data/CORUM_proteincomplexes.json', 'wb') as outfile:
            json.dump(corum_json, outfile)

        go_json = dict()
        for k,v in self.go.items():
            for label in v:
                go_json[label] = [ k ]

        with open('/Users/otvisitor/Documents/data/GO_macromolecularcomplexes.json', 'wb') as outfile:
            json.dump(go_json, outfile)


def main():
    pce_object = ProteinComplexExtractor()
    #pce_object.parse_ChEMBL()
    #pce_object.parse_CORUM()
    pce_object.parse_GO()
    pce_object.generate_bioentities_files()

if __name__ == "__main__":
    main()




