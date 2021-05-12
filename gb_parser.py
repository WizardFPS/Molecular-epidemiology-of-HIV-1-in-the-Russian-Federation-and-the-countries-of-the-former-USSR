import argparse
import csv
import os
import re
import sys
import time
from pathlib import Path

#input_file - genBank file

def parse_gb(input_file, min_size, max_size, suf):

    #File with countries/cities/other features and their abbreviations
    #Should be located at the same directory as script
    print(os.getcwd())
    COUNTRY_MAP_FILE = os.path.join(sys.path[0],"country_map.csv")
    FEATURE_MAP_FILE = os.path.join(sys.path[0],"feature_map.csv")
    CITY_MAP_FILE = os.path.join(sys.path[0],"rus_city_map.csv")

    #name of output file
    OUTPUT_FILE = '.'.join(input_file.split('.')[:-1]) + '.fasta'
    #length tresholds for sequences

    MIN_ORIGIN_SIZE = min_size
    MAX_ORIGIN_SIZE = max_size

    #qualifiers in GB entry which will added to name of fasta-sequence
    FIELD_NAMES = ["subtype", "sex","risk_factor","country","collection_date"]
    #fields from HIVDataBaseData
    COMMENT_NAMES = ["Infection country", "Participant sex", "Patient sex", "Mode of transmission", "Risk factor", "Subtype"]
    COUNTRIES = ["Russia", "Ukraine", "Belarus"]

    #country_map - &&&
    country_map = read_csv(COUNTRY_MAP_FILE)

    global countries_list
    print(country_map)
    '''
    countries_list = list(country_map)
    countries_list.sort()
    ind_nig = countries_list.index('niger')
    ind_nigeria = countries_list.index('nigeria')
    ind_gb = countries_list.index('guinea-bissau')
    ind_g = countries_list.index('guinea')
    ind_georg = countries_list.index('georgia')
    ind_usa = countries_list.index('usa')
    
    if (ind_nig and ind_nigeria):
        countries_list[ind_nigeria] = "niger"
        countries_list[ind_nig] = "nigeria"
    if (ind_gb and ind_g):
        countries_list[ind_g] = "guinea-bissau"
        countries_list[ind_gb] = "guinea"
    if (ind_nig and ind_nigeria):
        countries_list[ind_usa] = "georgia"
        countries_list[ind_georg] = "usa"
    '''
    #cl = list(country_map.keys())
    #print(cl)
    #cl.sort()
    #print(cl)

    city_map = read_csv(CITY_MAP_FILE)
    
    #feature_map = read_csv(FEATURE_MAP_FILE, strip_it=False)
    #city_map = read_csv(CITY_MAP_FILE)
    
    #RegExp that are used for parsing  INPUT_FILE
    accession = re.compile(r"^ACCESSION\s+([a-z_A-Z0-9]+)")
    definition = re.compile(r"^DEFINITION\s+([a-z_A-Z0-9]+)")
    #accession = re.compile(r"^LOCUS\s+([a-z_A-Z0-9]+)")
    #this searches for any qualifiers from FIELDS_NAMES regardless the feature_key
    features = re.compile(r"^\s+\/([a-z_]+)=\"([^\"]+)[\"]*")

    #this searches for field names from COMMENT feature_key
    features_hivdb_st = re.compile(r"[^\>]*##HIVDataBaseData-START##[^\>]*")
    features_hivdb_en =  re.compile(r"[^\>]*##HIVDataBaseData-END##")


    #subtype = re.compile(r"subtype[:=\s]*")
    
    #origin = re.compile(r"\s+\d+\s+([a-z ]+)$")
    origin = re.compile(r"\s+\d+\s+([atgcryswkmbdhvvn\s]+)$")
    end_line = re.compile(r"^//\s*$")
    #clean_from = re.compile(r"[:;.,/]")
    #clean_to = "_"
    clean_from = re.compile(r"[_\s,]")
    clean_to = "-"

    #Reg Exp for date
    #1994
    year0 =  re.compile(r"[0-9]{4}")
    year_between = re.compile(r"((collect(ed)*(ion)*)|between|during|date)([\s]+in)*([\s]+date[\s]+range)*[\s]+[0-9]{4}(\s)*(\-|and)(\s)*[0-9]{4}")
    #_1994_ 
    year1 = re.compile(r"[_\"\-\s,;/]+[0-9]{4}[_\"\-\s,;]+")
    #year11 = re.compile(r"[_\"\-]*[0-9]{4}[_\"\-]+")
    #isolated November 1994,isolated 1994, isolated on June 06,1994
    year2 = re.compile(r"(isolated|taken|sample(d)*|collect(ed)*(ion)*)[\s]*[a-zA-Z]*[\s]*[a-zA-Z]*[\s]*[0-9]*[\s,]*[0-9]{4}")
    #sample(d)* in 1994, sample 1994
    year3 = re.compile("sample(d)*[\s]+[a-zA-Z]*[\s]*[0-9]{4}")
    #collected/collection (in)* Nov* 1994
    year4 = re.compile("collect(ed)*(ion)*[\s,]+[a-zA-Z]*[\s,]*[a-zA-Z]*[\s,]*[0-9]{4}")
    #date November 1994
    year5 = re.compile(r"date[\s,:]*[A-Za-z]*[,\s]*[0-9]*[,\s]*[0-9]{4}")
    #02 Nov 1994
    year6 = re.compile(r"[\s,\"]+[0-9]{2}[\.\-\/\\\s]+[A-Za-z]{3}[\.\-\/\\\s]+[0-9]{4}[\s,;\.\"]+")
    #09-11-1994
    year7 =  re.compile(r"[\s,\"]+[0-9]{2}[\.\-\/\\\s]+[0-9]{2}[\.\-\/\\\s]+[0-9]{4}[\s,;\.\"]+")
    #1994-11-09
    year8 =  re.compile(r"[\s,\"]+[0-9]{4}[\.\-\/\\\s]+[0-9]{2}[\.\-\/\\\s]+[0-9]{2}[\s,;\.\"]+")
    #09-11-94
    year9 =  re.compile(r"[\s,\"]+[0-9]{2}[\.\-\/\\\s]+[0-9]{2}[\.\-\/\\\s]+[0-9]{2}[\s,;\.\"]+")

    def check_year(dict_data, qual='', stri=''):
        if qual != '':
            stri = dict_data[qual]
        m = year_between.search(stri)
        if m:
            year = re.search(r"[0-9]{4}(\s)*\-(\s)*[0-9]{4}",m.group())
            if year:
                dict_data["collection_date"] = year.group()
                return
        for reg_year in [year1,year2,year3,year4,year5,year6,year7,year8]:
            m = reg_year.search(stri)
            if m:
                year_st = year0.search(m.group())
                if year_st:
                    year = int(year_st.group())
                    if year > 1940 and year < 2018:
                        dict_data["collection_date"] = str(year)
                        return
        '''
        m = year9.search(stri)
        if m:
            year = int(m.group().strip('\"').strip(' ').strip(';').strip('\.')[-2:])
            if year>60:
                dict_data["collection_date"] ='19'+str(year)
                return
            else:
                if year<=18:
                    dict_data["collection_date"] ='20'+str(year)
                    return
        '''
        
        for month in list_months:
            m1 = re.search(month+r"[\s,.]+[0-9]{2}",stri)
            m2 = re.search(month+r"[\s,.]+[0-9]{4}",stri)
            if m1:
                year = int(m1.group()[-2:])
                if year>60:
                    dict_data["collection_date"] ='19'+str(year)
                    return
                else:
                    if year<=18:
                        dict_data["collection_date"] ='20'+str(year)
                        return
            if m2:
                year_st = m2.group()[-4:]
                year = int(year_st)
                if year > 1940 and year < 2018:
                    dict_data["collection_date"] = str(year)
                    return

        '''
        m = re.search(r"[\s,;\-\"]+[0-9]{4}[\s,;\-\"]+",stri)
        if m:
            year_st = year0.search(m.group()).group()
            if year_st:
                year = int(year_st)
                if year > 1940 and year < 2018:
                    dict_data["collection_date"] = str(year)
                    return
        '''
        return


    #checks if the string dict_data[qual] contains any key from dict_par 
    def check_keys(dict_par, par, dict_data, qual='', st = ''):
        for key in dict_par.keys():
            if qual == '':
                if dict_par == dict_risk:
                    r_exp =re.compile(r"[Ff]actor[)(\s:,;_\-\"]+{0}[)(\s,;_\"\n]+".format(key))
                else:
                    r_exp =re.compile(r"[)(\s:,;_\-\"]+{0}[)(\s,;_\"\n]+".format(key))
            else:
                r_exp =re.compile(r"[)(\s:,;_\-\"]+{0}[)(\s,;_\"\n]+".format(key))
            
            #r_exp =re.compile(r"[^\]>+{0}[\s,;_\"\n]+".format(key))
            if qual !='':
                 st = dict_data[qual].lower()

            m = r_exp.search(st.lower())
            if m:
                dict_data[par] = dict_par[key]
                return
        return

    def check_subtype(dict_data, qual='', stri=''):

        subtype = re.compile(r"(inter)*[Ss]ub[-]*type[:=\s]+([a-zA-Z0-9_\/\-]+\s*)+(\([a-zA-Z0-9_\/\-]+\)\(*[a-zA-Z0-9_\/\-]*\)*)*")
        genotype = re.compile(r"[Gg]enotype[:=\s]+([a-zA-Z0-9_\/\-]+\s*)+(\([a-zA-Z0-9_\/\-]+\)\(*[a-zA-Z0-9_\/\-]*\)*)*")
        subcluster = re.compile(r"[Ss]ubcluster[:=\s]+([a-zA-Z0-9_\/\-]+\s*)+(\([a-zA-Z0-9_\/\-]+\)\(*[a-zA-Z0-9_\/\-]*\)*)*")
        clade = re.compile(r"[Cc]lade[:=\s]+([a-zA-Z0-9_\/\-]+\s*)+(\([a-zA-Z0-9_\/\-]+\)\(*[a-zA-Z0-9_\/\-]*\)*)*")
        mixed = re.compile(r"[Mm]ixed\s*[Ss]ubtype[s]*[:=\s]+[pP]rotease=[a-zA-Z_0-9]+[;,\s]+[Rr][Tt]=[a-zA-Z_0-9]+")

        if qual != '':
            stri =dict_data[qual]

        for reg in [mixed,subtype, genotype, subcluster,clade]:
            m = reg.search(stri)
            if m != None:
                break

        st=None
        if m:
            all = m.group(0)
            pr = re.search(r"[pP]rotease=[a-zA-Z_0-9]+[;,\s]+", all)
            rt = re.search(r"[rR][tT]=[a-zA-Z0-9_]+", all)
            if rt and pr:
                st = re.search(r"=[a-zA-Z_0-9]+",pr.group()).group().strip('=') + re.search(r"=[a-zA-Z_0-9]+",rt.group()).group().strip('=')
                dict_data["subtype"] = re.sub(clean_from, clean_to, st)
            else:
                a = re.search(r"\([a-zA-Z0-9_\/\-]+\)\(*[a-zA-Z0-9_\/\-]*\)*", all)
                if a:
                    st = re.sub("\)",".", re.sub("\(","",a.group(0).strip(')').strip('(')))
                    dict_data["subtype"] = re.sub(clean_from, clean_to, st)
                    return
                else:
                    b = re.findall(r"[a-zA-Z0-9_\/\-]+",all)
                    if b[1] == 'rec' or b[1] == 'Rec' or b[1] == 'Recomb' or b[1] == 'recombinant' or b[1] == 'HIV-1':
                        if len(b) ==3:
                            st = re.sub(clean_from, clean_to, b[2])
                            return
                        else:
                            st='rec'
                    else:
                        if b[1] == 'CRF' or b[1] == 'crf' or b[1] == 'Crf':
                            st= re.sub(clean_from, clean_to, b[1]+b[2])
                        else:
                            st = re.sub(clean_from, clean_to, b[1])
                    dict_data["subtype"] = st
                return
        else:
            a = re.match(r'["]+[a-zA-Z_]+[0-9]*\n["]*', stri)
            if a:
                st = a.group().strip('\n').strip('"')
                dict_data["subtype"] = re.sub(clean_from, clean_to, st)
                return
            else:
                #m = re.search(r'[NPOM]:[A-DF-HJK][0-9]*',stri)
                m = re.search(r'[M]:[A-DF-HJK][0-9]*',stri)
                if m:
                    st=re.sub(clean_from, clean_to, m.group())
                    dict_data["subtype"] = re.sub(clean_from, clean_to, st[2:])
                    return
                else:
                    m = re.search(r'[UC]RF\s*[0-9A-Za-z_\-]*', stri)
                    if m:
                         st = re.sub(clean_from, clean_to, m.group()) 
                         dict_data["subtype"] = st
                         return
        return


    def check_country(stri, test_features):
        for key in dict_country:
            if re.match('\"'+key.lower(),stri.lower()):
                test_features["country"]=dict_country[key]
        return


    #Dictionaries with values for sex, transmission route
    dict_sex = {'female': 'F',
                'sex[\s:]+[fF]': 'F',
                'woman': 'F',
                'girl':'F',
                'boy':'M',
                'male': 'M',
                'sex[\s:]+[Mm]': 'M',
                'man': 'M',
                'transexual':'T'
                }
    
    dict_risk = {     #intravenous drug use
                      'idu':'IDU',
                     'iv[\s]+drug[\s]+use':'IDU',
                     'idus':'IDU',
                     'ivdu':'IDU',
                     'ivdus':'IDU',
                     'people[\s]+who[\s]+inject[\s]+drugs':'IDU',
                     'drug[\s]+use[r]*[s]*':'IDU',
                     'pi':'IDU',
                     #male sex with male
                     'm[ae]n[\s\-]+who[\s\-]+have[\s\-]+sex[\s\-]+with[\s\-]+m[ae]n':'MSM',
                     'm[ae]n[\s\-]+hav[e]*(ing)*[\s\-]+sex[\s\-]+with[\s\-]+m[ae]n':'MSM',
                     'male[\s-]+sex[\s\-]+with[\s\-]+male':'MSM',
                     'msm': 'MSM',
                     'sm': 'MSM',
                     'homosexual': 'MSM',
                     'homosexually-infected':'MSM',
                     'sg':'MSM',
                     'male[\s\-]+to[\s\-]+male[\s\-]+transmission':'MSM',
                     #heterosexual
                     'heterosexual':'HET',
                     'hetero':'HET',
                     'hsx':'HET',
                     'sh':'HET',
                     'het':'HET',
                     'hts':'HET',
                     'mswm':'SEX',
                     'sexual':'SEX',
                     #sex worker
                     'sw':'SEX',
                     #bisexual
                     'sb':'SEX',
                     'bisexual':'SEX',
                      #sexually undescribed
                     'su':'SEX',
                     'male[\s\-]+to[\s\-]+female[\s\-]+transmission':'HET',
                     'female[\s\-]+to[\s\-]+male[\s\-]+transmission':'HET',
                     'vaginal[\s]+intercourse':'HET', 
                     'insertive[\s]+partner':'HET',
                     'receptive[\s]+partner': 'HET',
                     'heterosexually':'HET',
                     'heterosexually-infected':'HET',
                     'heterosexsual':'HET',
                     'idu-heterosexual':'IDU-HET',
                     'ver':'VER',
                     'vertical[\s]+transmission':'VER',
                     'infected[\s]+infant':'VER',
                     'infected[\s]+child':'VER',
                     'mother[\s\-]+to[\s\-]+child[\s]+transmission':'VER',
                     'child[\s]*born[\s]*to[\s]*HIV[\-1\s]*infected[\s]*mother':'VER',
                     'born[\s]*to[\s]*HIV[\-1\s]*infected[\s]*woman':'VER',
                     'HIV[\-1\s]*infected[\s]*child':'VER',
                     'mtct':'VER',
                     'perinatally':'VER',
                     'mb':'VER',
                     #nosocomial
                     'no':'NKM',
                     'nosocomial':'NKM',
                     #blood transfusion
                     'pb':'BTR',
                     'factor[\s]+blood':'BTR',
                     '(blood|transfusion)[\s]+(donor|recipient)':'BTR',
                     'blood[\s]+transfusion':'BTR',
                     'par[ea]nteral[l]*[y]*':'BTR',
                     'medical[\s]+manipulation':'MED',
                     #experimental
                     'ex':'EX',
                     'sex':'SEX'
                     #'not reported':'NR'
                     }
    dict_country = {
        'ukrainian':'UKR',
        'russia': 'RUS',
        'ru':'RUS',
        'rus':'RUS',
        'rus':'RUS',
        'lat': 'LVA',
        'est':'EST',
        'lit':'LTU',
        'ukranian':'UKR',
        'ua':'UKR',
        'blr':'BLR',
        'br':'BRA',
        'ar':'ARG',
        'gf':'FGU',
        'col':'COL'}
    
    list_months = ['january','february','march','april','may','june','july','august','september','sep','october','november','december']

    #Parsing INPUT_FILE
    
    #all data is stored into list test
    tests = []
    #entries with no collection date
    tests_nodate = []
    
    #accession number
    test_accession = ""
    test_definition= ""
    #all features which are written like /***=""
    test_features = {}
    test_hivdatabase = ''
    #'origin' field
    test_origin = ""
    #try:
    list_ac=[]


    #qualifier of key_feature we are looking for
    current_key = ''
    print(input_file)
    input_dir = '\\'.join(input_file.split('\\')[:-2])+'\\'
    print(input_dir+'debuging\\'+suf+'\\note.txt')
    
    
    Path(input_dir+'debuging\\'+suf).mkdir(parents=True, exist_ok=True)
    
    
    
    f_note =open(input_dir+'debuging\\'+suf+'\\note.txt', 'w')
    f_organism =open(input_dir+'debuging\\'+suf+'\\organism.txt', 'w')
    f_isols =open(input_dir+'debuging\\'+suf+'\\isols.txt', 'w')
    f_host =open(input_dir+'debuging\\'+suf+'\\host.txt', 'w')
    f_isol =open(input_dir+'debuging\\'+suf+'\\isol.txt', 'w')
    f_strain =open(input_dir+'debuging\\'+suf+'\\strain.txt', 'w')
    f_hivdb =open(input_dir+'debuging\\'+suf+'\\hivdb.txt', 'w')

    #flag for source feature key
    s=0
    #number of string
    num=0
    #flag for qualifier inside source
    j=0
    #flag for HIVDataBaseData section
    g=0

    #number of entry
    num_entry = 0
    #definition flag
    d=0
    for line in open(input_file, "r"):
        if num % 10000 == 0:
            print(num)
        #print(num)
        num+=1
        m = definition.match(line)
        #finds DEFINITION field
        if(m):
            test_definition = line.lower().strip('\n')
            d=1
            continue

        #finds ACCESSION field using RegExp
        #print(line)
        m = accession.match(line)
        if(m):
            test_accession = m.group(1)
            list_ac.append(test_accession)
            num_entry+=1
            d=0
        if d==1:
            test_definition = test_definition + ' ' + line.lower().strip('\n')

        
        
        #HIVDataBaseData

        ##finds ##HIVDataBaseDataStart
        m = features_hivdb_st.match(line)
        if(m):
            #flag indicates that we are in HIVDataBaseData section
            g=1
            a= re.sub(r"\s+"," ",m.group())
            test_hivdatabase = test_hivdatabase+ a.strip('\n') + ' '
            continue

        ##finds ##HIVDataBaseDataEnd
        m = features_hivdb_en.match(line)
        if(m):
            g=0
            a= re.sub(r"^\s+"," ",m.group())
            test_hivdatabase = test_hivdatabase+ a + ' '
            continue
        if g==1:
            test_hivdatabase =  re.sub(r"^\s+"," ",test_hivdatabase)+ line.strip('\n') + ' '
            continue

        #qualifiers in /source feature key

        #check if /source has started or has ended
        if re.match(r"^\s{5}[a-zA-z_]+[\s0-9<>\.]+", line) or re.match(r"^[a-zA-z]+", line):
            s=0
        if re.match(r"^\s+source[0-9\.\s]+", line):
            s=1
        
        if s==1:
            #finds all qualifiers using RegExp
            m = features.match(line)

            if(m):
                o = re.search("\"",m.group(0)).start()
                test_features[m.group(1)] = m.group(0)[o:].strip('\n')
                current_key = m.group(1)
                j=1
                continue

            if j==1:
                test_features[current_key] = test_features[current_key] + ' '+ re.sub("^\s+"," ",line.strip('\n'))
                

                #print(m.group(1))
                #print(m.group(2))

        #ORIGIN field - dna sequence is written in this field
        m = origin.match(line)
            
        if(m):
            test_origin += re.sub("\s+","",m.group(1))

        if(end_line.search(line)):
            
            s=0
            g=0
            j=0

            if num_entry %1000 == 0:
                print('\t'+str(num_entry))
            #print('\t'+str(num_entry))
            t1 = time.time()

            #check HIVDataBaseData for sex, subtype, risk factor, country
            if test_hivdatabase != '':
                if "country" not in test_features:
                    map_country(test_hivdatabase, country_map)
                if "sex" not in test_features:
                    check_keys(dict_sex, "sex", test_features, st=test_hivdatabase)
                if "risk_factor" not in test_features:
                    check_keys(dict_risk, "risk_factor", test_features, st=test_hivdatabase)
                if "subtype" not in test_features:
                    check_subtype(test_features, stri=test_hivdatabase)
                if "collection_date" not in test_features:
                    check_year(test_features, stri=test_hivdatabase)
            t2 = time.time()
            #print('t2-t1='+str(t2-t1))
            '''
            if "Participant sex" in test_comments:
                if test_comments["Participant sex"] != 'no date':
                    test_features['sex']=dict_sex[test_comments["Participant sex"].lower()]

            if "Patient sex" in test_comments:
                if test_comments["Patient sex"] != 'no date':
                    test_features['sex']=dict_sex[test_comments["Patient sex"].lower()]

            if "Risk factor" in test_comments:
                if test_comments["Risk factor"] != 'no date':
                    test_features['risk_factor'] = dict_risk[test_comments["Risk factor"].lower()]

            if "Subtype" in test_comments:
                #test_features["subtype"]=
                if test_comments["Subtype"] != 'no date':
                    test_features["subtype"]= test_comments["Subtype"]
                    print(test_features["subtype"])
            '''

            if ("collection_date" or "collected_by") in test_features: 
                year = year0.search(test_features["collection_date"])
                if year:
                    test_features["collection_date"] =  year.group()
                else:
                    year = year9.search(test_features["collection_date"])
                    if year:
                        year = int(year.group().strip("\"")[-2:])
                        if year>60:
                            test_features["collection_date"] ='19'+str(year)
                        else:
                            test_features["collection_date"] ='20'+str(year)
            t3 = time.time()
            #print('t3-t2='+str(t3-t2))
            if "country" in test_features:
                country = map_country(test_features["country"].strip('\"'), country_map, fl=0)
                if country == "RUS":
                    city = map_city(test_features["country"], city_map)
                    if city !=None:
                        print(test_features["country"], city)
                        test_features["country"] = country + ':'+city
                        test_features["city"] = city
                    else:
                        test_features["country"] = country
                else:
                    if country == 'GBR':
                        print(test_accession, test_features["country"])
                    test_features["country"] = country
            if "sex" in test_features:
                check_keys(dict_sex, "sex", test_features, qual="sex")
            if "genotype" in test_features:
                if re.search(r'[a-zA-Z0-9_\-]+',test_features["genotype"]):
                    test_features["subtype"] = re.sub(clean_from, clean_to,test_features["genotype"].strip('"'))
            if "subtype" in test_features:
                test_features["subtype"] = re.sub(clean_from, clean_to, test_features["subtype"].strip('"'))

            t4 = time.time()
            #print('t4-t3='+str(t4-t3))
            qualifiers = ["note", "isolation_source", "host", "strain", "isolate", "organism"]
            for qual in qualifiers:
                if qual in test_features:
                    if "sex" not in test_features:
                        check_keys(dict_sex, "sex", test_features, qual=qual)
                    if "risk_factor" not in test_features:
                        check_keys(dict_risk, "risk_factor", test_features, qual=qual)
                    if "subtype" not in test_features:
                        check_subtype(test_features, qual=qual)
                    if "collection_date" not in test_features:
                        check_year(test_features,qual=qual)
            t5 = time.time()
            #print('t5-t4='+str(t5-t4))
            '''
            for key in ["sex", "risk_factor"]:
                if key not in test_features:
                    if "isolate" in test_features:
                        if key == "sex":
                            check_keys(dict_sex, "sex", test_features, "isolate")
                        if key == "risk_factor":
                            check_keys(dict_risk, "risk_factor", test_features, "isolate")
                    if key not in test_features:
                        if "strain" in test_features:
                            if key == "sex":
                                check_keys(dict_sex, "sex", test_features, "strain")
                            if key == "risk_factor":
                                check_keys(dict_risk, "risk_factor", test_features, "strain")
            '''
            #fetch country from "strain" field


            t52 =  time.time()
            #print('t52-t51='+str(t52-t51))
            for qual in ["note", "isolate", "strain"]:
                if "country" not in test_features:
                    if qual in test_features:
                         country = map_country(test_features[qual],country_map)
                         if country:
                            test_features["country"] = country
                            if country == 'GBR':
                                print(test_accession)
            if "country" not in test_features:
                country = map_country(test_definition, country_map)
                if country:
                    test_features["country"] = country



            t51 =  time.time()
            #print('t51-t5='+str(t51-t5))
            #fetch country from "isolate" filed
            if "country" not in test_features:
                if "isolate" in test_features:
                    check_country(test_features["isolate"], test_features)
            if "country" not in test_features:
                if "strain" in test_features:
                    check_country(test_features["strain"], test_features)

            t6 = time.time()
            #if num_entry %1000 == 0:
                #print('t6-t52='+str(t6-t52))
            
            if "sex" not in test_features:
                if "risk_factor" in test_features:
                    if test_features["risk_factor"] == 'MSM':
                        test_features["sex"]='M'
            t7 = time.time()
            #print('t7-t6='+str(t7-t6))
            #if "subtype" not in test_features:
            #    if "isolate" in test_features:
            #        m = re.match(r'[A-DF-HJK][0-9]*', test_features['isolate'].split('_')[-1])
            #        if m:
            #            test_features["subtype"] = m.group(0)
            if ("country" in test_features) and (test_features["country"]=="RUS"):
                if "city" not in test_features:
                    for qual in ["note","isolation_source"]:
                        if qual in test_features:
                            city = map_city(test_features[qual].strip('\"'), city_map)
                            if city !=None:
                                print(city)
                                test_features["country"] = test_features["country"] + ':'+city
                                test_features["city"] = city
                                break


            if "host" in test_features:
                f_host.write(test_accession+ ' '+test_features["host"]+'\n')
            if "organism" in test_features:
                f_organism.write(test_accession+ ' '+test_features["organism"]+'\n')
            if "isolate" in test_features:
                f_isol.write(test_accession+ ' '+test_features["isolate"]+'\n')
            if "isolation_source" in test_features:
                f_isols.write(test_accession+ ' '+test_features["isolation_source"]+'\n')
            if "note" in test_features:
                f_note.write(test_accession+ ' '+test_features["note"]+'\n')
            if "strain" in test_features:
                f_strain.write(test_accession+ ' '+test_features["strain"]+'\n')
            if test_hivdatabase!='':
                f_hivdb.write(test_accession+ ' '+test_hivdatabase+'\n')



            #for k in ["strain", "organism"]:
                #if k in test_features: 
                    #test_features[k] = map_feature(test_features[k], feature_map)


            if ("collection_date" or "collected_by") not in test_features:
                tests_nodate.append({
                    "accession" : test_accession, 
                    "features" : test_features,
                    "origin" : test_origin
                    })
                    
            #else:
                #if test_features["collection_date"] == '2016' or test_features["collection_date"] == '2017':
                    #print(test_accession)
                    
            else:
                tests.append({
                    "accession" : test_accession, 
                    "features" : test_features,
                    "origin" : test_origin
                    })

            test_accession = ""
            test_definition = ""
            test_features = {}
            test_hivdatabase = ''
            test_origin = ""
    #except:
    #   print('Error: can\'t read genbank file')

    f_note.close()
    f_organism.close()
    f_isols.close()
    f_host.close()
    f_isol.close()
    f_strain.close()
    f_hivdb.close()

    # Write results into OUTPUT_FILE 

    #print(tests_nodate)
    
    for ent in tests_nodate:
        tests.append(ent)

    out = open(OUTPUT_FILE, "w+")

    data_table_f = ".".join(input_file.split(".")[:-1]) + ".csv"
    data_table = open(data_table_f, "w")
    data_table.write('Index,GB_ID,subtype,sex,risk_factor,country,year\n')
    for test in tests:
        aaccession, f, origin = test["accession"], test["features"], test["origin"]
        #print(aaccession, f, origin)

        if not (MIN_ORIGIN_SIZE < len(origin) < MAX_ORIGIN_SIZE):
            #print(aaccession, len(origin))
            continue

        output = ">%s" % aaccession
        output_table = aaccession
        for field_name in FIELD_NAMES:
            if field_name not in f:
                output += "_%s" % "NA"
                output_table += ",%s" % "NA"
                continue
            output += "_%s" % f[field_name]
            output_table += ",%s" % f[field_name]
        output_table =output.strip('>')+','+output_table+ "\n"
        output += "\n"
        


        #out.write(re.sub(clean_from, clean_to, output))
        out.write(output)
        out.write(origin + "\n")

        data_table.write(output_table)
        
    out.close()
    data_table.close()
    return list_ac


def read_csv(file_name, strip_it=True):
    if not os.path.exists(file_name):
        return {}

    def strip(value):
        return strip_it and value.strip() or value

    with open(file_name) as csvfile:
        reader = csv.DictReader(csvfile,
                                delimiter=",", 
                                fieldnames=["base", "new"])

        result = {}
        for row in reader:
            result[strip(row["base"].lower())] = strip(row["new"])
        return result


def map_country(country, country_map, fl=0):
    country_lower =  country.lower()
    t_en = time.time()
    #for k in countries_list:
    for k in country_map.keys():
        t_st = time.time()
        if country_lower == k:
            return country_map[k]
        else:
            
            #if country_lower.startswith(k) or search or search1 or search2:
            if  country_lower.startswith(k):
                return country_map[k]
            else:
                search = re.search(r"[a-z]*[,\s.\"]+"+k.lower(),country_lower)
                search1 = re.search(r"country[\s]+"+k.lower(),country_lower)
                if search or search1:
                #print('+t_en-t_st='+str(t_en-t_st)+' '+country)
                    return country_map[k]

    t_en = time.time()
    #print('t_en-t_st='+str(t_en-t_st))
    if fl ==1:
        return country
    else:
        return None

    #return country

def map_feature(feature, feature_map):
    feature_lower = feature.lower()

    for k, v in feature_map.items():
        return re.sub(k, v, feature, flags=re.I)

    return feature

#replaces the name of city by its coordinates

def map_city(city, city_map):
    city_lower =  city.lower()
    for k, v in city_map.items():
        reg_exp = r'[\s;,]+'+k.lower()+'[\s;,\"]+'
        if re.search(reg_exp, city_lower):
            #print(re.search(reg_exp, city_lower).group(),v)
            return v

    return None
#l1 = parse_gb('D://MY_FILES//DATA//Lukashev//HIV//FSU_041218.gb')
#l1.sort()
#l2 = parse_gb('D://MY_FILES//DATA//Lukashev//HIV//hiv-db.genbank')
#l2.sort()
#set_1, set_2 = set(l1), set(l2)
#print(len(list(set_1 & set_2)))

#file_name = sys.argv[1]
#suffix = sys.argv[2]

#file_name =  'D://MY_FILES//DATA//Lukashev//HIV//gb//HIVdbfsu250419.genbank'
#suffix = 'fsu'

#print(file_name)
#parse_gb(file_name,suffix)
#parse_gb('D://MY_FILES//DATA//Lukashev//HIV//gb//Europe_part1.genbank')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required = True)
    parser.add_argument("-min", "--min_length", type=int,
                        help="Minimal length of sequence.\
                        Sequences shorter than min length will not be included in the final dataset",
                         required = True)
    parser.add_argument("-max", "--max_length", type=int,
                        help="Maximal length of sequence. \
                        Sequences longer than max length will not be included in the final dataset",
                         required = True)
    parser.add_argument("-suf", "--suffix", type=str,
                        help="Name of folder where certain string for debugging will be saved",
                         required = True)
    args = parser.parse_args()

    #if not len(sys.argv) == 9:
    #    print("Please, use \"python parser_gb.py --help\"")
    #else:
    parse_gb(args.input_file, args.min_length, args.max_length, args.suffix)