"""
Functions for Simulations for Pangenomic Analyses

Is was based on the LabBCES XYLMAN project:
https://github.com/labbces/CAZymes-MAN-XYL

"""

# create connection object to mysql/mariadb database using sqlalchemy
from sqlalchemy.sql.sqltypes import Enum, String


def connectDB(password=None):
    import sqlalchemy
    from sqlalchemy import create_engine
    #connect to postgresql database
    db_user = 'genomicdatabase'
    db_password = password
    db_name = 'genomedb'
    db_host = 'localhost'
    db_port =  3306
    db_url = 'mariadb+pymysql://{}:{}@{}:{}/{}?use_unicode=1&charset=utf8'.format(db_user,db_password,db_host,db_port,db_name)
    engine = create_engine(db_url,encoding='utf-8')
    return engine


def createDB(password=None):
    """Create database schema.


    """
    engine = connectDB(password)
    from sqlalchemy.ext.declarative import declarative_base
    from sqlalchemy.orm import relationship
    from sqlalchemy import Column, Integer, String, ForeignKeyConstraint, PrimaryKeyConstraint, Text, Date, ForeignKey, UniqueConstraint

    Base = declarative_base()

    class Taxonomy(Base):
        """Add docstrings

        """
        __tablename__ = 'Taxonomy'
        #id = Column(Integer, primary_key=True)
        ParentTaxID = Column(Integer, nullable=True)
        TaxID = Column(Integer, primary_key=True, autoincrement=False)
        #TaxIDRank = Column(Integer)
        RankName = Column(String(255))
        TaxName = Column(String(255))
        __table_args__ = (
            PrimaryKeyConstraint(TaxID),
            {'mariadb_engine':'InnoDB'},
        )

    class Genome(Base):
        """Add docstrings

        """
        __tablename__ = 'Genomes'
        __table_args__ = (
            ForeignKeyConstraint(['TaxID'], ['Taxonomy.TaxID']),
            {'mariadb_engine':'InnoDB'},
                    )
        AssemblyAccession = Column(String(100), primary_key=True)
        TaxID = Column(Integer)
        urlBase = Column(String(255))

    class GenomeFile(Base):
        """Add docstrings
        
        """
        __tablename__ = 'GenomeFiles'
        ID=Column(Integer, primary_key=True, autoincrement=True)
        AssemblyAccession = Column(String(100))
        FileSource = Column(String(255))
        FileType = Column(String(100))
        FileName = Column(String(255))
        __table_args__ = (
            ForeignKeyConstraint(['AssemblyAccession'], ['Genomes.AssemblyAccession']),
            UniqueConstraint(FileType, AssemblyAccession),
            {'mariadb_engine':'InnoDB'},
                    )

    class GenomeFileDownloaded(Base):
        """Add docstrings
        
        """
        __tablename__ = 'GenomeFileDownloaded'
        ID=Column(Integer, primary_key=True, autoincrement=True)
        GenomeFileID = Column(Integer)
        Action=Column(String(255))
        ActionDate = Column(Date)
        __table_args__ = (
            ForeignKeyConstraint(['GenomeFileID'], ['GenomeFiles.ID']),
            UniqueConstraint(GenomeFileID, Action, ActionDate),
            {'mariadb_engine':'InnoDB'},
        )

    class Protein2GenomeFile(Base):
        """Add docstrings
        
        """
        __tablename__ = 'Proteins2GenomeFile'
        ID=Column(Integer, primary_key=True, autoincrement=True)
        GenomeFileID = Column(Integer, index=True)
        ProteinID=Column(String(255))
        __table_args__ = (
            ForeignKeyConstraint(['GenomeFileID'], ['GenomeFiles.ID']),
            ForeignKeyConstraint(['ProteinID'], ['ProteinSequences.ProteinID']),
            UniqueConstraint(GenomeFileID,ProteinID),
            {'mariadb_engine':'InnoDB'},
        )

    class ProteinSequence(Base):
        """Add docstrings

        """
        __tablename__ = 'ProteinSequences'
        ProteinID = Column(String(255), primary_key=True)
        Database =Column(String(50))
        Sequence = Column(Text)
        HashSequence = Column(String(255))

    Base.metadata.create_all(engine)


#Drop all tables from DB
def dropDB(password=None):
    engine = connectDB(password)
    from sqlalchemy.ext.automap import automap_base
    Base = automap_base()
    Base.prepare(engine, reflect=True)
    Base.metadata.drop_all(engine)


def computeMD5Sumfile(pathFile=None):
    import hashlib
    with open(pathFile, "rb") as f:
        file_hash = hashlib.md5()
        chunk = f.read(8192)
        while chunk:
            file_hash.update(chunk)
            chunk = f.read(8192)
    return file_hash.hexdigest()


def getMD5sumFromFile(md5PathFile=None, target=None):
    with open(md5PathFile, "r") as fmd5:
        for line in fmd5:
            line = line.rstrip()
            fields = line.split()
            if fields[1].replace('./','') == target:
                md5sum = fields[0]
                break
    return md5sum


def getTaxInfo(taxID=None,password=None):
    engine = connectDB(password)
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy.sql import text
    from sqlalchemy import select, update, bindparam
    from sqlalchemy.ext.automap import automap_base
 
    Base = automap_base()
    Base.prepare(engine, reflect=True)
    Session = sessionmaker(bind=engine)
    session = Session()

    taxonomicLineage = text("""with recursive ancestors as ( 
         select * from Taxonomy where TaxID=:taxID 
             union 
             select f.*  from Taxonomy as f, ancestors as a 
             where  f.TaxID = a.ParentTaxID 
             )
             select a.TaxName, a.RankName, a.TaxID 
             from ancestors as a, Taxonomy as b 
             where a.ParentTaxID = b.TaxID""")
    resultsGetTaxonomicLineage=session.execute(taxonomicLineage, {'taxID':taxID})
    rows=resultsGetTaxonomicLineage.fetchall()
    lineage={}
    if rows:
        targetGroup='Do not know - Problem with taxonomy'
        name=''
        for row in rows:
            # print(f'{row[0]} {row[1]} {row[2]}')
            lineage[row[1]]=row[0]
            if row[1]=='superkingdom' and row[0]=='Bacteria':
                targetGroup=row[0]
            elif row[1]=='superkingdom' and row[0]=='Archaea':
                targetGroup=row[0]
            elif row[1]=='kingdom' and row[0]=='Fungi':
                targetGroup=row[0]
            elif row[2] == taxID:
                name=row[0]
        lineage['targetGroup']=targetGroup
        lineage['name']=name
    return lineage


def downloadGenomeFiles(password=None, dirPath=None, fileType=None):
    engine = connectDB(password)
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy import select, update, bindparam
    from sqlalchemy.ext.automap import automap_base
    import sys
    from os import path, mkdir, makedirs,remove
    import time
    import urllib.request
    import datetime
    import platform

    Base = automap_base()
    Base.prepare(engine, reflect=True)
    Session = sessionmaker(bind=engine)
    session = Session()

    Genomes = Base.classes.Genomes
    GenomeFiles = Base.classes.GenomeFiles
    GenomeFileDownloaded = Base.classes.GenomeFileDownloaded

    if path.exists(dirPath):
        print(f'Directory \'{dirPath}\' exists. OK.')	
    else:
        print(f'Directory \'{dirPath}\' does not exists. Creating it.')
        mkdir(dirPath)
    #Get the list of genomes that have not been downloaded yet, 50 at the time (recursive)
    getListOfGenomes2Download=select([Genomes.AssemblyAccession,Genomes.urlBase,GenomeFiles.FileName,GenomeFiles.ID]).join(GenomeFiles,Genomes.AssemblyAccession==GenomeFiles.AssemblyAccession).join(GenomeFileDownloaded,GenomeFiles.ID==GenomeFileDownloaded.GenomeFileID,isouter=True).where(GenomeFileDownloaded.ID==None).where(GenomeFiles.FileType=='Protein sequence').limit(100)
    print(getListOfGenomes2Download)
    resultGetListOfGenomes2Download=session.execute(getListOfGenomes2Download)
    rows=resultGetListOfGenomes2Download.fetchall()
    counterFiles2Download=0
    if(rows):
        for row in rows:
            subDirs=row[1].replace('https://ftp.ncbi.nlm.nih.gov/genomes/all/','')
            if platform.system() == 'Windows':
                subDirs=subDirs.replace('/','\\')
            elif platform.system() == 'Linux':
                subDirs=subDirs.replace('\\','/')
            urlFile=row[1]+'/'+row[2]
            completePath=path.join(dirPath,subDirs)
            completePathFile=path.join(completePath,row[2])
            completePathMD5File=path.join(completePath,'md5checksums.txt')
            urlMD5File=row[1]+'/md5checksums.txt'
            print(f'{urlFile}')
            if path.exists(completePath):
                print(f'Directory \'{completePath}\' exists. OK.')	
            else:
                print(f'Directory \'{completePath}\' does not exists. Creating it.')
                makedirs(completePath)
            if path.exists(completePathFile):
                print(f'Genome file {row[2]} exists', file=sys.stderr)
                md5sumLocalFile=computeMD5Sumfile(completePathFile)
                if path.exists(completePathMD5File):
                    md5sumRemoteFile=getMD5sumFromFile(completePathMD5File,row[2])
                else:
                    urllib.request.urlretrieve(urlMD5File, completePathMD5File)
                    md5sumRemoteFile=getMD5sumFromFile(completePathMD5File,row[2])
            else:
                counterFiles2Download+=1
                if counterFiles2Download % 10 == 0:
                    time.sleep(1)
                    session.commit()
                print(f'Genome  file {row[2]}  does not exist,  start downloading', file=sys.stdout)
                urllib.request.urlretrieve(urlFile, completePathFile)
                urllib.request.urlretrieve(urlMD5File, completePathMD5File)
                md5sumLocalFile=computeMD5Sumfile(completePathFile)
                md5sumRemoteFile=getMD5sumFromFile(completePathMD5File,row[2])
            if md5sumLocalFile != md5sumRemoteFile:
                print(f'MD5 checksum of file {row[2]} does not match. Deleting file.', file=sys.stderr)
                remove(completePathFile)
            else:
                print("Download complete", file=sys.stderr)
                getGenomeFileDownloaded=select([GenomeFileDownloaded.GenomeFileID]).where(GenomeFileDownloaded.Action=='Downloaded').where(GenomeFileDownloaded.GenomeFileID==row[3])
                resultsGetGenomeFileDownloaded=session.execute(getGenomeFileDownloaded)
                if resultsGetGenomeFileDownloaded.fetchone() is None:
                    print(f'Genome file {row[2]} downloaded. Inserting action \'Downloaded\' and date into DB.', file=sys.stdout)
                    dateToday=datetime.date.today()
                    session.add(GenomeFileDownloaded(GenomeFileID=row[3], Action='Downloaded', ActionDate=dateToday))
                else:
                    print(f'Genome file {row[2]} downloaded. Action \'Downloaded\' already in DB.', file=sys.stdout)
    else:
        sys.exit('No more files to process')

    session.commit()
    session.close()
    downloadGenomeFiles(password=password, dirPath=dirPath, fileType=fileType)


def updateProteinSequences(password=None,apiKey=None):
    """Populate the ProteinSequences table.
    
    Updates by getting the protein sequences from the Genbank/Uniprot
    """
    engine = connectDB(password)
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy import select, update, bindparam
    from sqlalchemy.ext.automap import automap_base
    import sys

    Base = automap_base()
    Base.prepare(engine, reflect=True)
    Session = sessionmaker(bind=engine)
    session = Session()

    ProteinSequences = Base.classes.ProteinSequences

    updateStmt = (
        update(ProteinSequences).
        where(ProteinSequences.ProteinID == bindparam('proteinID')).
        where(ProteinSequences.Database == bindparam('DB')).
        values(Sequence=bindparam('sequence'))
    )

    getProteinSequences=select([ProteinSequences.ProteinID,ProteinSequences.Database]).where(ProteinSequences.Sequence==None).limit(100)
    resultsGetProteinSequences=session.execute(getProteinSequences)
    rows=resultsGetProteinSequences.fetchall()
    proteinIDsGenbank=[]
    proteinIDsUniprot=[]
    if rows:
        for row in rows:
            if row[1]=='genbank':
                proteinIDsGenbank.append(row[0])
                #sequence = getProteinSequenceFromGenbank(proteinID=row[0],apiKey=apiKey)
            elif row[1]=='uniprot':
                # sequence = getProteinSequenceFromUniprot(row[0])
                proteinIDsUniprot.append(row[0])
            else:
                print(f'Database {row[1]} not supported yet',file=sys.stderr)
        if len(proteinIDsGenbank)>0:
            # print(f'Updating {len(proteinIDsGenbank)} {proteinIDsGenbank} protein sequences from genbank',file=sys.stderr)	
            seqsGenbank=getProteinSequenceFromGenbank(proteinIDs=proteinIDsGenbank,apiKey=apiKey)
            # print(seqsGenbank,file=sys.stderr)
            if seqsGenbank:
                session.execute(updateStmt, seqsGenbank)
        if len(proteinIDsUniprot)>0:
            # print(f'Updating {len(proteinIDsUniprot)} protein sequences from uniprot',file=sys.stderr)
            seqsUniprot=getProteinSequenceFromUniprot(proteinIDs=proteinIDsUniprot)
            if seqsUniprot:
                session.execute(updateStmt, seqsUniprot)
                # session.commit()
    else:
        sys.exit('No more proteins to process')
        
    session.commit()
    session.close()
    updateProteinSequences(password=password, apiKey=apiKey)
    

def getProteinSequenceFromUniprot(proteinIDs=None):
    """Gets seqeunces from UniProt.

    """
    import sys
    from io import StringIO
    from Bio import SeqIO
    import requests
    import time
    
    seqsList=[]
    counter=0
    baseUrl="http://www.uniprot.org/uniprot/"

    for id in proteinIDs:
        counter+=1
        if counter%100==0:
            time.sleep(3)
        currentUrl=baseUrl+id+".fasta?version=*"
        #print(currentUrl)
        response = requests.post(currentUrl)
        cData=''.join(response.text)
        if cData:
            fastaIO=StringIO(cData)
            seqObj=SeqIO.parse(fastaIO,'fasta')
            seqsDict={}
            seqsDict['proteinID']=id
            seqsDict['DB']='uniprot'
            # print(seqsDict)
            seqsDict['sequence']=str(next(seqObj).seq)
            # print(seqsDict)
            # print(seqObj)
            seqsList.append(seqsDict)
            fastaIO.close()
        else:
            print(f'No sequence found for {id}',file=sys.stderr)
    return seqsList
    

def getProteinSequenceFromGenbank(proteinIDs, apiKey=None):
    """Gets sequences from Genbank.
    
    """
    import sys
    from io import StringIO
    from Bio import Entrez, SeqIO
    import re
    import time
    
    seqsList=[]
    if(apiKey):
        Entrez.email = "diego.riano@cena.usp.br"
        Entrez.api_key = apiKey
        searchRes = Entrez.read(Entrez.epost("protein", id=",".join(proteinIDs)))
        webenv = searchRes["WebEnv"]
        query_key = searchRes["QueryKey"]
        
        try:
            fastaIO = StringIO(Entrez.efetch(db="protein", rettype="fasta", retmode="text", webenv=webenv, query_key=query_key, api_key=apiKey).read())
        except:
            print(f'Error while getting sequences from Genbank for {proteinIDs}',file=sys.stderr)
            return seqsList
        
        seqsObj=SeqIO.parse(fastaIO,'fasta')
        for seq in seqsObj:
            # print(f'XXX:{seq.id}',file=sys.stderr)
            match=re.search(r'^[a-z]*\|+([A-Z0-9.]*)(\|.+)?$',seq.id,re.IGNORECASE)
            if match:
                # print('match')
                #Dealing with weird cases, where the accession retrieved is different. Some specific cases.
                if match.group(1) == '3PPS':
                    proteinID='333361328'
                elif match.group(1) == '4K3A':
                    proteinID='550545166'
                else:
                    #this deals with th enormal case, when the id of the retrieve seq should match wit the ID stored in the DB
                    proteinID=match.group(1)
            else:
                proteinID=seq.id
            seqsDict={}
            sequence=str(seq.seq)
            #if proteinID not in seqsDict:
            #    seqsDict[proteinID]={}
            seqsDict['proteinID']=proteinID
            seqsDict['DB']='genbank'
            seqsDict['sequence']=sequence
            seqsList.append(seqsDict)
            #updateProteinSequences(password=None,apiKey=apiKey,proteinID=proteinID,sequence=sequence)
        fastaIO.close()
        time.sleep(3)
        return seqsList
    else:
        print(f'No API key provided for Entrez. Please provide one in the command line.',file=sys.stderr)


def populateGenomes(url,password=None,updateNCBITaxDB=False,typeOrg='euk'):
    """Add docstrings
    
    """
    engine = connectDB(password)

    from sqlalchemy.orm import sessionmaker
    from sqlalchemy import select
    from sqlalchemy.ext.automap import automap_base
    from ete3 import NCBITaxa
    import urllib.request
    import ftplib
    import time
    from os.path import exists, getmtime
    import sys

    ncbi = NCBITaxa()
    
    if updateNCBITaxDB:#TODO. This is not working. Check it.
        ncbi.update_taxonomy_database()
    
    if typeOrg == 'euk':
        assemblyAccessionIndex=8
    else:
        assemblyAccessionIndex=18

    counter=0
    #create a session
    Base = automap_base()
    Base.prepare(engine, reflect=True)
    Session = sessionmaker(bind=engine)
    session = Session()
    
    Genome = Base.classes.Genomes
    Taxonomy = Base.classes.Taxonomy
    GenomeFile = Base.classes.GenomeFiles

    #get the data from the NCBI genomes and add to tables
    if exists('genomes.txt'):
        print(f'Genome info file from NCBI exists,  checking age', file=sys.stderr)
        mtime=getmtime('genomes.txt')
        now=time.time()
        if (now - mtime) > 604800: #Download the file again if the file in disk is older than 7 days (60*60*24*7)
            print(f'Genome info file from NCBI is older than 7 days, downloading', file=sys.stderr)
            urllib.request.urlretrieve(url, 'genomes.txt')
            print("Download complete", file=sys.stderr)
        else:
            print(f'Genome info file from NCBI is younger than 7 days, not downloading and processing as it is', file=sys.stderr)
    else:
        print(f'Genome info file from NCBI does not exist,  start downloading', file=sys.stderr)
        urllib.request.urlretrieve(url, 'genomes.txt')
        print("Download complete", file=sys.stderr)

    with open('genomes.txt', mode='r', encoding="utf8") as genomeFile, open('populateErrors.log','w') as errorsLog:
        for line in genomeFile:
            #line = line.decode('utf-8')
            line = line.rstrip()
            fields = line.split('\t')
            if line.startswith('#'):
                continue
            else:
                #If the genome info was already inserted skip it
                checkGenome=select([Genome]).where(Genome.AssemblyAccession==fields[assemblyAccessionIndex])
                resultCheckGenome = session.execute(checkGenome)
                if resultCheckGenome.fetchone():
                    # print(f'Already saw {fields[assemblyAccessionIndex]}')
                    continue
                else:
                    True
                    # print(f'Checking {fields[assemblyAccessionIndex]}')

                #Commit to the DB every 100 lines of the genome file
                if counter % 100 == 0:
                    session.commit()

                #If there is no taxonomy info for the genome, control the error and continue
                try:
                    lineage = ncbi.get_lineage(int(fields[1]))
                except:
                    errorsLog.write(f'Error: Missing TaxID from NCBI DB: {fields[1]} for assembly: {fields[assemblyAccessionIndex]}\n')
                    continue
                #Only processes genomes with a taxonomy ID in fungi, archaea or bacteria
                targetGroups={4751,2157,2,10239}
                if targetGroups.intersection(set(lineage)):
                    counter+=1
                    print(f'\tAdding Genome Info for Assembly {fields[assemblyAccessionIndex]}')
                    if(lineage[-1] != fields[1]):
                        errorsLog.write(f'Warning: TaxID from genome file: {fields[1]} was translated to: {lineage[-1]} for assembly: {fields[assemblyAccessionIndex]}\n')
                    for index, i in enumerate(lineage):
                        checkTaxID=select([Taxonomy]).where(Taxonomy.TaxID==i)
                        resultCheckTaxID = session.execute(checkTaxID)
                        if i == 1:
                            #Check whether the taxonomy info is already in the DB. If not, add it
                            if resultCheckTaxID.fetchone() is None:
                                taxid2name = ncbi.get_taxid_translator([i])
                                rank = ncbi.get_rank([i])
                                # print(f'{index}\t{i}\t{lineage[index-1]}\t{lineage[index]}\t{taxid2name[i]}\t{rank[i]}')
                                #The root node is added to the DB, as it does not have parent and NA will be inserted for ParentTaxID
                                session.add(Taxonomy(TaxID=i, RankName=rank[i], TaxName=taxid2name[i])) 
                        else:
                            #Check whether the taxonomy info is already in the DB. If not, add it
                            if resultCheckTaxID.fetchone() is None:
                                taxid2name = ncbi.get_taxid_translator([i])
                                rank = ncbi.get_rank([i])
                                # print(f'{index}\t{i}\t{lineage[index-1]}\t{lineage[index]}\t{taxid2name[i]}\t{rank[i]}')
                                session.add(Taxonomy(ParentTaxID=lineage[index-1], TaxID=i, RankName=rank[i], TaxName=taxid2name[i]))
                    checkGenome=select([Genome]).where(Genome.AssemblyAccession==fields[assemblyAccessionIndex])
                    resultCheckGenome = session.execute(checkGenome)
                    #Check which files are available for the genome.
                    if resultCheckGenome.fetchone() is None:
                        acc,ver=fields[assemblyAccessionIndex].split('.')
                        if len(acc) == 13:
                            FTP_USER = "anonymous"
                            FTP_PASS = ""
                            FTP_HOST = "ftp.ncbi.nlm.nih.gov"
                            try:
                                ftp = ftplib.FTP(FTP_HOST, FTP_USER, FTP_PASS)
                                ftp.cwd(f'genomes/all/{acc[0:3]}/{acc[4:7]}/{acc[7:10]}/{acc[10:13]}/')
                                for asm in ftp.nlst():
                                    if asm.startswith(fields[assemblyAccessionIndex]):
                                        url=f'https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[0:3]}/{acc[4:7]}/{acc[7:10]}/{acc[10:13]}/{asm}'
                                        # print(f'{acc}\t{ver}\t{asm}\t{url}')
                                        #Sometimes the NCBI can update TaxIDs, i.e., translate the taxID. NCBITaxa deals with this, but to avoid foreigkey errors
                                        #it is better to insert into the DB what NCBITaxa retrieved and not the TaxID in the genome file
                                        session.add(Genome(AssemblyAccession=fields[assemblyAccessionIndex], TaxID=lineage[-1], urlBase=url))
                                        ftp.cwd(asm)
                                        for file in ftp.nlst():
                                            if file.endswith(asm + '_genomic.fna.gz'):
                                                session.add(GenomeFile(AssemblyAccession=fields[assemblyAccessionIndex], FileType='Genome sequence', FileName=file, FileSource='NCBI'))
                                            elif file.endswith(asm + '_protein.faa.gz'):
                                                session.add(GenomeFile(AssemblyAccession=fields[assemblyAccessionIndex], FileType='Protein sequence', FileName=file, FileSource='NCBI'))	
                                            elif file.endswith(asm + '_genomic.gff.gz'):
                                                session.add(GenomeFile(AssemblyAccession=fields[assemblyAccessionIndex], FileType='Genome annotation', FileName=file, FileSource='NCBI'))
                                            elif file.endswith(asm + '_translated_cds.faa.gz'):
                                                session.add(GenomeFile(AssemblyAccession=fields[assemblyAccessionIndex], FileType='Protein sequence alter', FileName=file, FileSource='NCBI'))
                                            # print(file)
                                ftp.quit()
                                time.sleep(3)#Wait to avoid being banned by NCBI
                            except:
                                errorsLog.write(f'Error: Could not download files for assembly: {fields[assemblyAccessionIndex]}\n')
                else:
                    print(f'\t{fields[assemblyAccessionIndex]} with TaxID: {fields[1]} is not in the target groups')
    #commit the changes
    session.commit()
    #close the session
    session.close()