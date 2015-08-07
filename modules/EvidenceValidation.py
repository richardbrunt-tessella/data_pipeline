import os
import sys
import re
import gzip
import smtplib
import time
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.MIMEBase import MIMEBase
from email import Encoders
import logging
from StringIO import StringIO
import json
from json import JSONDecoder
from json import JSONEncoder
from datetime import datetime
from sqlalchemy import and_, or_, table, column, select, update, insert
from common import Actions
from common.PGAdapter import *
import cttv.model.core as cttv
import cttv.model.flatten as flat
#import cttv.model.flatten as flat
from settings import Config
import hashlib
BLOCKSIZE = 65536


__author__ = 'gautierk'

logger = logging.getLogger(__name__)

# figlet -c "Validation Passed"
messagePassed='''
__     __    _ _     _       _   _               ____                        _
\ \   / /_ _| (_) __| | __ _| |_(_) ___  _ __   |  _ \ __ _ ___ ___  ___  __| |
 \ \ / / _` | | |/ _` |/ _` | __| |/ _ \| '_ \  | |_) / _` / __/ __|/ _ \/ _` |
  \ V / (_| | | | (_| | (_| | |_| | (_) | | | | |  __/ (_| \__ \__ \  __/ (_| |
   \_/ \__,_|_|_|\__,_|\__,_|\__|_|\___/|_| |_| |_|   \__,_|___/___/\___|\__,_|
   
'''

messageFailed='''
   __     __    _ _     _       _   _               _____     _ _          _
   \ \   / /_ _| (_) __| | __ _| |_(_) ___  _ __   |  ___|_ _(_) | ___  __| |
    \ \ / / _` | | |/ _` |/ _` | __| |/ _ \| '_ \  | |_ / _` | | |/ _ \/ _` |
     \ V / (_| | | | (_| | (_| | |_| | (_) | | | | |  _| (_| | | |  __/ (_| |
      \_/ \__,_|_|_|\__,_|\__,_|\__|_|\___/|_| |_| |_|  \__,_|_|_|\___|\__,_|

'''

cttv_data_providers_e_mails = {
 "cttv001" : [ 'gautierk@targetvalidation.org', 'mmaguire@ebi.ac.uk', 'samiulh@targetvalidation.org', 'andreap@targetvalidation.org' ],
 "cttv006" : [ 'fabregat@ebi.ac.uk' ],
 "cttv007" : [ 'kl1@sanger.ac.uk' ],
 "cttv008" : [ 'mpaulam@ebi.ac.uk', 'patricia@ebi.ac.uk' ],
 "cttv009" : [ 'gautierk@targetvalidation.org', 'mmaguire@ebi.ac.uk'], #[ 'cleroy@ebi.ac.uk' ],
 "cttv010" : [ 'mkeays@ebi.ac.uk' ],
 "cttv011" : [ 'eddturner@ebi.ac.uk' ],
 "cttv012" : [ 'fjlopez@ebi.ac.uk', 'garys@ebi.ac.uk' ],
 "cttv025" : [ 'kafkas@ebi.ac.uk', 'ftalo@ebi.ac.uk' ] 
}

efo_current = {}
efo_obsolete = {}
ensembl_current = {}

class EvidenceValidationActions(Actions):
    CHECKFILES='checkfiles'
    VALIDATE='validate'

class EvidenceValidationFileChecker():

    def __init__(self, adapter):
        self.adapter = adapter
        self.session = adapter.session
        #formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        #self.buffer = StringIO()
        #streamhandler = logging.StreamHandler(self.buffer)
        #streamhandler.setFormatter(formatter)
        #memoryhandler = logging.handlers.MemoryHandler(1024*10, logging.DEBUG, streamhandler)        
        #LOGGER = logging.getLogger('cttv.model.core')
        #LOGGER.setLevel(logging.ERROR)
        #LOGGER.addHandler(memoryhandler)
        #LOGGER = logging.getLogger('cttv.model.evidence.core')
        #LOGGER.setLevel(logging.ERROR)
        #LOGGER.addHandler(memoryhandler)        
        #LOGGER = logging.getLogger('cttv.model.evidence.association_score')
        #LOGGER.setLevel(logging.ERROR)
        #LOGGER.addHandler(memoryhandler)

    def startCapture(self, newLogLevel = None):
        """ Start capturing log output to a string buffer.
        
        http://docs.python.org/release/2.6/library/logging.html
        
        @param newLogLevel: Optionally change the global logging level, e.g. logging.DEBUG
        """
        self.buffer = StringIO()
        self.logHandler = logging.StreamHandler(self.buffer)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        formatter = logging.Formatter("%(asctime)s - %(message)s")
        self.logHandler.setFormatter(formatter)
        
        #print >> self.buffer, "Log output"
        
        #for module in [ 'cttv.model.core', 'cttv.model.evidence.core', 'cttv.model.evidence.association_score' ]:
        rootLogger = logging.getLogger()
     
        if newLogLevel:
            self.oldLogLevel = rootLogger.getEffectiveLevel()
            rootLogger.setLevel(newLogLevel)
        else:
            self.oldLogLevel = None
        
        rootLogger.addHandler(self.logHandler)    
            
    def stopCapture(self):
        """ Stop capturing log output.
        
        @return: Collected log output as string        
        """
                                
        # Remove our handler
        #for module in [ 'cttv.model.core', 'cttv.model.evidence.core', 'cttv.model.evidence.association_score' ]:
        rootLogger = logging.getLogger()        

        # Restore logging level (if any)
        if self.oldLogLevel:
            rootLogger.setLevel(self.oldLogLevel)

        rootLogger.removeHandler(self.logHandler)
        
        self.logHandler.flush()
        self.buffer.flush()
        
        return self.buffer.getvalue()
    
    def send_email(self, bSend, provider_id, filename, bValidated, nb_records, errors, when, extra_text, logfile):
        me = "support@targetvalidation.org"
        you = ",".join(cttv_data_providers_e_mails[provider_id])
        status = "passed"
        if not bValidated:
            status = "failed"
        # Create message container - the correct MIME type is multipart/alternative.
        #msg = MIMEMultipart('alternative')
        msg = MIMEMultipart()
        msg['Subject'] = "CTTV: Validation of submitted file {0} {1}".format(filename, status)
        msg['From'] = me
        msg['To'] = you
        rcpt = cttv_data_providers_e_mails[provider_id]
        if provider_id != 'cttv001':
            rcpt.extend(cttv_data_providers_e_mails['cttv001'])
            msg['Cc'] = ",".join(cttv_data_providers_e_mails['cttv001'])

        text = "This is an automated message generated by the CTTV Core Platform Pipeline on {0}\n".format(when)

        if bValidated:
            text += messagePassed
            text += "Congratulations :)\n"
        else:
            text += messageFailed
            text += "See details in the attachment {0}\n\n".format(os.path.basename(logfile))
        text += "JSON schema version:\t1.2.1\n"
        text += "Number of evidence strings:\t{0}\n".format(nb_records)
        for key in errors:
            text += "Number of {0}:\t{1}\n".format(key, errors[key])
        text += "\n"
        text += extra_text
        text += "\nYours\nThe CTTV Core Platform Team"
        #text += signature
        print text
        
        if bSend:
            # Record the MIME types of both parts - text/plain and text/html.
            part1 = MIMEText(text, 'plain')
            
            # Attach parts into message container.
            # According to RFC 2046, the last part of a multipart message, in this case
            # the HTML message, is best and preferred.
            msg.attach(part1)
            
            if not bValidated:
                part2 = MIMEBase('application', "octet-stream")
                part2.set_payload( open(logfile,"rb").read() )
                Encoders.encode_base64(part2)
                part2.add_header('Content-Disposition', 'attachment; filename="{0}"'.format(os.path.basename(logfile)))
                msg.attach(part2)        
            
            # Send the message via local SMTP server.
            mail = smtplib.SMTP('smtp.office365.com', 587)

            mail.ehlo()

            mail.starttls()

            mail.login(me, 'P@ssword')
            mail.sendmail(me, rcpt, msg.as_string())
            mail.quit()
        return 0;
    
    def load_Ensembl(self):
        logging.info("Loading Ensembl {0} asssembly genes".format(Config.EVIDENCEVALIDATION_ENSEMBL_ASSEMBLY))
        for row in self.session.query(EnsemblGeneInfo).filter_by(assembly_name = Config.EVIDENCEVALIDATION_ENSEMBL_ASSEMBLY).all():
            #print "%s %s"%(row.ensembl_gene_id, row.external_name)
            ensembl_current[row.ensembl_gene_id] = row.external_name
            
    def load_efo(self):
        # Change this in favor of paths
        logging.info("Loading EFO current and obsolete terms")
        efo_pattern = "http://www.ebi.ac.uk/efo/EFO_%"
        orphanet_pattern = "http://www.orpha.net/ORDO/Orphanet_%"
        hpo_pattern = "http://purl.obolibrary.org/obo/HP_%"
        mp_pattern = "http://purl.obolibrary.org/obo/MP_%"
        
        for row in self.session.query(EFONames).filter(
                or_(
                    EFONames.uri.like(efo_pattern), 
                    EFONames.uri.like(mp_pattern), 
                    EFONames.uri.like(orphanet_pattern), 
                    EFONames.uri.like(hpo_pattern)
                    )
                ):
            #print row.uri
            efo_current[row.uri] = row.label
        for row in self.session.query(EFOObsoleteClass):
            #print "obsolete %s"%(row.uri)
            efo_obsolete[row.uri] = row.reason
            
    def check_all(self):
    
        self.load_Ensembl();
        self.load_efo();
        #return;
        
        for dirname, dirnames, filenames in os.walk(Config.EVIDENCEVALIDATION_FTP_SUBMISSION_PATH):
            for subdirname in dirnames:
                #cttv_match = re.match("^(cttv[0-9]{3})$", subdirname)
                cttv_match = re.match("^(cttv001)$", subdirname)
                if cttv_match:
                    provider_id = cttv_match.groups()[0]
                    cttv_dir = os.path.join(dirname, subdirname)
                    print(cttv_dir)
                    for cttv_dirname, cttv_dirs, filenames in os.walk(os.path.join(cttv_dir, "upload/submissions")):
                        for filename in filenames:
                            if filename.endswith(('.json.gz')) and filename == 'cttv_external_disgenet-29-07-2015.json.gz': #and filename == 'cttv009-14-07-2015.json.gz': #
                                cttv_file = os.path.join(cttv_dirname, filename)
                                last_modified = os.path.getmtime(cttv_file)
                                july = time.strptime("01 Jul 2015", "%d %b %Y")
                                julyseconds = time.mktime(july)
                                if( last_modified - julyseconds ) > 0:
                                    m = re.match("^(.+).json.gz$", filename)
                                    logfile = os.path.join(cttv_dirname, m.groups()[0] + "_log.txt")
                                    print(cttv_file)
                                    md5_hash = self.check_gzipfile(cttv_file)
                                    self.validate_gzipfile(cttv_file, filename, provider_id, md5_hash, logfile = logfile)

        self.session.commit()
        
    def validate_gzipfile(self, file_on_disk, filename, provider_id, md5_hash, logfile = None):
        '''
        check if the file was already processed
        '''
        bValidate = False
        bGivingUp = False
        rowToUpdate = None
        count = self.session.query(EvidenceValidation.filename).filter_by(filename=file_on_disk).count()
        logging.info('Was the file parsed already? %i'%(count))
        if count == 0:
            bValidate = True
        else:
            for row in self.session.query(EvidenceValidation).filter_by(filename=file_on_disk):
                if row.md5 == md5_hash:
                    logging.info('%s == %s'% (row.md5, md5_hash))
                    logging.info('%s file already recorded. Won\'t parse'%file_on_disk)
                    return;
                else:
                    logging.info('%s != %s'% (row.md5, md5_hash))
                    bValidate = True
                    rowToUpdate = row
                    break;
        logging.info('bValidate %r'% (bValidate))
        # Check EFO overrepresentation
        # Check target overrepresentation
        diseases = {}
        top_diseases = []
        targets = {}
        top_targets = []
        obsolete_diseases = {}
        invalid_diseases = {}
        invalid_ensembl_ids = {}
        
        if bValidate == True:
            logging.info('Starting validation of %s'% (file_on_disk))
            fh = gzip.GzipFile(file_on_disk, "r")
            #lfh = gzip.open(logfile, 'wb', compresslevel=5)
            lfh = open(logfile, 'wb')
            cc = 0
            lc = 0
            nb_errors = 0
            nb_duplicates = 0
            nb_efo_invalid = 0
            nb_efo_obsolete = 0
            nb_ensembl_invalid = 0
            hexdigest_map = {}
            for line in fh:
                #logging.info(line)
                python_raw = json.loads(line)
                # now validate
                obj = None
                validation_result = 0
                validation_failed = False
                
                if ('label' in python_raw  or 'type' in python_raw) and 'validated_against_schema_version' in python_raw and python_raw['validated_against_schema_version'] == "1.2.1":
                    if 'label' in python_raw:
                        python_raw['type'] = python_raw.pop('label', None)
                    data_type = python_raw['type']
                    #logging.info('type %s'%data_type)
                    if data_type in ['genetic_association', 'rna_expression', 'genetic_literature', 'affected_pathway', 'somatic_mutation', 'known_drug', 'literature', 'animal_model']:
                        if data_type == 'genetic_association':
                            obj = cttv.Genetics.fromMap(python_raw)
                        elif data_type == 'rna_expression':
                            obj = cttv.Expression.fromMap(python_raw)
                        elif data_type in ['genetic_literature', 'affected_pathway', 'somatic_mutation']:
                            obj = cttv.Literature_Curated.fromMap(python_raw)
                        elif data_type == 'known_drug':
                            obj = cttv.Drug.fromMap(python_raw)
                            #logging.info(obj.evidence.association_score.__class__.__name__)
                            #logging.info(obj.evidence.target2drug.association_score.__class__.__name__)
                            #logging.info(obj.evidence.drug2clinic.association_score.__class__.__name__)
                        elif data_type == 'literature':
                            obj = cttv.Literature_Mining.fromMap(python_raw)
                        elif data_type == 'animal_model':
                            obj = cttv.Animal_Models.fromMap(python_raw)

                        if obj.target.id:
                            for id in obj.target.id: 
                                if id in targets:
                                    targets[id] +=1
                                else:
                                    targets[id] = 1
                                if not id in top_targets:
                                    if len(top_targets) < Config.EVIDENCEVALIDATION_NB_TOP_TARGETS:
                                        top_targets.append(id)
                                    else:
                                        # map,reduce
                                        for n in range(0,len(top_targets)):
                                            if targets[top_targets[n]] < targets[id]:
                                                top_targets[n] = id;
                                                break;
                                        
                        if obj.disease.id:
                            for id in obj.disease.id:
                                if id in diseases:
                                    diseases[id] +=1
                                else:
                                    diseases[id] =1
                                if not id in top_diseases:
                                    if len(top_diseases) < Config.EVIDENCEVALIDATION_NB_TOP_DISEASES:
                                        top_diseases.append(id)
                                    else:
                                        # map,reduce
                                        for n in range(0,len(top_diseases)):
                                            if diseases[top_diseases[n]] < diseases[id]:
                                                top_diseases[n] = id;
                                                break;
                                        
                        if not bGivingUp:  
                            self.startCapture(logging.ERROR)
                        uniq_elements = obj.unique_association_fields
                        uniq_elements_flat = flat.DatatStructureFlattener(uniq_elements)
                        uniq_elements_flat_hexdig = uniq_elements_flat.get_hexdigest()
                        if not uniq_elements_flat_hexdig in hexdigest_map:
                            hexdigest_map[uniq_elements_flat_hexdig] = [ lc+1 ]
                        else:
                            hexdigest_map[uniq_elements_flat_hexdig].append(lc+1)                          
                            logger.error("Line {0}: Duplicated unique_association_fields on lines {1}".format(lc+1, ",".join(map(lambda x: "%i"%x,  hexdigest_map[uniq_elements_flat_hexdig]))))
                            nb_duplicates = nb_duplicates + 1
                            validation_failed = True
  
                        validation_result = obj.validate(logger)
                        nb_errors = nb_errors + validation_result
                        
                        '''
                        Check EFO
                        '''
                        if obj.disease.id:
                            for disease_id in obj.disease.id:
                                if disease_id not in efo_current:
                                    logger.error("Line {0}: Invalid disease term detected {1}. Please provide the correct EFO disease term".format(lc+1, disease_id))
                                    if disease_id not in invalid_diseases:
                                        invalid_diseases[disease_id] = 1
                                    else:
                                        invalid_diseases[disease_id] += 1
                                    nb_efo_invalid +=1
                                if disease_id in efo_obsolete:
                                    logger.error("Line {0}: Obsolete disease term detected {1} ('{2}'): {3}".format(lc+1, disease_id, efo_current[disease_id], efo_obsolete[disease_id]))
                                    if disease_id not in obsolete_diseases:
                                        obsolete_diseases[disease_id] = 1
                                    else:
                                        obsolete_diseases[disease_id] += 1
                                    nb_efo_obsolete +=1
                                
                        '''
                        Check Ensembl ID, UniProt ID (TODO)
                        '''
                        if obj.target.id:
                            for id in obj.target.id:
                                # http://identifiers.org/ensembl/ENSG00000178573
                                m = re.match('http://identifiers.org/ensembl/(ENSG\d+)', id)
                                if m:
                                    ensembl_id = m.groups()[0].rstrip("\s")
                                    if not ensembl_id in ensembl_current:
                                        logger.error("Line {0}: Invalid Ensembl gene detected {1}. Please provide the correct identifier for assembly {2}".format(lc+1, ensembl_id, Config.EVIDENCEVALIDATION_ENSEMBL_ASSEMBLY))
                                        if not ensembl_id in invalid_ensembl_ids:
                                            invalid_ensembl_ids[ensembl_id] = 1
                                        else:
                                            invalid_ensembl_ids[ensembl_id] += 1
                                        nb_ensembl_invalid +=1
                                    
                                
                        if not bGivingUp:  
                            logs = self.stopCapture()
                    else:
                        if not bGivingUp:
                            self.startCapture(logging.ERROR)
                            logger.error("Line {0}: '{1}' is not a valid 1.2.1 evidence string type".format(lc+1, data_type))
                            logs = self.stopCapture()
                        nb_errors += 1
                        validation_failed = True
                elif not 'validated_against_schema_version' in python_raw or ('validated_against_schema_version' in python_raw and python_raw['validated_against_schema_version'] != "1.2.1"):
                    if not bGivingUp:
                        self.startCapture(logging.ERROR)
                        logger.error("Line {0}: Not a valid 1.2.1 evidence string - please check the 'validated_against_schema_version' mandatory attribute".format(lc+1))
                        logs = self.stopCapture()
                    nb_errors += 1
                    validation_failed = True
                else:
                    if not bGivingUp:
                        self.startCapture(logging.ERROR)
                        logger.error("Line {0}: Not a valid 1.2.1 evidence string - please add the mandatory 'type' attribute".format(lc+1))
                        logs = self.stopCapture()
                    nb_errors += 1
                    validation_failed = True
                    
                if (validation_failed or  validation_result > 0) and not bGivingUp:
                    if obj:
                        lfh.write("line {0} - {1}".format(lc+1, json.dumps(obj.unique_association_fields)))
                    else:
                        lfh.write("line {0} ".format(lc+1))
                    lfh.write(logs)
                    if nb_errors > Config.EVIDENCEVALIDATION_MAX_NB_ERRORS_REPORTED or nb_duplicates > Config.EVIDENCEVALIDATION_MAX_NB_ERRORS_REPORTED:
                        lfh.write("Too many errors: giving up.\n")
                        bGivingUp = True
                    
                lc += 1
                cc += len(line)
            logging.info('nb line parsed %i (size %i)'% (lc, cc))
            fh.close()
            lfh.close()
            
            # write top diseases / top targets
            text = ""
            if top_diseases:
                text +="Top %i diseases:\n"%(Config.EVIDENCEVALIDATION_NB_TOP_DISEASES)
                for n in range(0,len(top_diseases)):
                    if top_diseases[n] in efo_current:
                        text +="\t-{0}:\t{1} ({2:.1f}%) {3}\n".format(top_diseases[n], diseases[top_diseases[n]], diseases[top_diseases[n]]*100/lc, efo_current[top_diseases[n]])
                    else:
                        text +="\t-{0}:\t{1} ({2:.1f}%)\n".format(top_diseases[n], diseases[top_diseases[n]], diseases[top_diseases[n]]*100/lc)
                text +="\n"
            if top_targets:
                text +="Top %i targets:\n"%(Config.EVIDENCEVALIDATION_NB_TOP_TARGETS)
                for n in range(0,len(top_targets)):
                    text +="\t-{0}:\t{1} ({2:.1f}%)\n".format(top_targets[n], targets[top_targets[n]], targets[top_targets[n]]*100/lc)
                text +="\n"

            # report invalid/obsolete EFO term
            if nb_efo_invalid > 0:
                text +="%i distinct invalid EFO terms found in %i (%.1f%s) of the records:\n"%(len(invalid_diseases), nb_efo_invalid, nb_efo_invalid*100/lc, '%' )
                for disease_id in invalid_diseases:
                    text += "\t%s\t(reported %i times)\n"%(disease_id, invalid_diseases[disease_id])
                text +="\n"
            if nb_efo_invalid > 0:
                text +="%i obsolete EFO terms found in %i (%.1f%s) of the records:\n"%(len(obsolete_diseases), nb_efo_obsolete, nb_efo_obsolete*100/lc, '%' )
                for disease_id in obsolete_diseases:
                    text += "\t%s\t(reported %i times)\t%s\n"%(disease_id, obsolete_diseases[disease_id], efo_obsolete[disease_id].replace("\n", " "))
                text +="\n"

            # report invalid Ensembl genes
            if nb_ensembl_invalid > 0:
                text +="%i distinct invalid Ensembl identifiers found in %i (%.1f%s) of the records:\n"%(len(invalid_ensembl_ids), nb_ensembl_invalid, nb_ensembl_invalid*100/lc, '%' )
                for ensembl_id in invalid_ensembl_ids:
                    text += "\t%s\t(reported %i times)\n"%(ensembl_id, invalid_ensembl_ids[ensembl_id])
                text +="\n"
                
            now = datetime.utcnow()

            if count == 0:
                # insert
                f = EvidenceValidation(
                    provider_id = provider_id, 
                    filename = file_on_disk,
                    md5 = md5_hash,
                    date_created = now,
                    date_modified = now,
                    date_validated = now,
                    nb_submission = 1,
                    nb_records = lc,
                    nb_errors = nb_errors,
                    nb_duplicates = nb_duplicates,
                    successfully_validated = (nb_errors == 0 and nb_duplicates == 0)
                )
                self.session.add(f)
                logging.info('inserted %s file in the validation table'%file_on_disk)
            else:
                # update database
                rowToUpdate.md5 = md5_hash
                rowToUpdate.nb_records = lc
                rowToUpdate.nb_errors = nb_errors
                rowToUpdate.nb_duplicates = nb_duplicates
                rowToUpdate.date_modified = now
                rowToUpdate.date_validated = now
                rowToUpdate.successfully_validated = (nb_errors == 0 and nb_duplicates == 0)
                self.session.add(rowToUpdate)
    
            self.send_email(
                True, 
                provider_id, 
                filename, 
                (nb_errors == 0 and nb_duplicates == 0 and nb_efo_invalid == 0 and nb_efo_obsolete == 0 and nb_ensembl_invalid == 0), 
                lc, 
                { 'JSON errors': nb_errors, 'duplicates': nb_duplicates, 'invalid EFO terms': nb_efo_invalid, 'obsolete EFO terms': nb_efo_obsolete, 'invalid Ensembl ids': nb_ensembl_invalid }, 
                now, 
                text, 
                logfile
                )
                                
    
    def check_gzipfile(self, filename):
    
        hasher = hashlib.md5()
        with gzip.open(filename,'rb') as afile:
        #with open(filename, 'rb') as afile:
            buf = afile.read(BLOCKSIZE)
            while len(buf) > 0:
                hasher.update(buf)
                buf = afile.read(BLOCKSIZE)
        md5_hash = hasher.hexdigest() 
        print(md5_hash)
        afile.close()
        return md5_hash

                