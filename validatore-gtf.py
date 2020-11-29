# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 20:28:14 2020

@author: Michele Milesi 844682
"""

import re

#Riceve come parametro la stringa che contiene il nome del file
#e legge riga per riga
def read_file(file_name):
    with open(file_name, 'r') as input_file:
        file_rows = input_file.readlines()
    return file_rows

#Prende come argomento un dizionario, la chiave e il valore da aggiungere
#alla chiave
def update_dict(dictionary, key, value):
    value_list = dictionary.get(key, [])
    value_list.append(value)
    dictionary[key] = value_list
    return

#Prende come argomento una lista di stringhe
#rimuove i caratteri di spazio alla fine di ogni stringa,
#rimuove tutti i commenti e scarta le righe vuote
#inoltre associa ad ogni riga la sua posizione (numero di riga) all'interno del file
#restituisce la nuova lista di tuple (numero_riga, riga)
def format_rows(string_list):
    format_string_list = []
    for count in range(0, len(string_list)):
        row = string_list[count]
        format_row = re.search('((("[^"]*")|[^#])*)#?', row.rstrip()).group(1)
        if format_row != '':
            format_string_list.append(((count+1), format_row))
    return format_string_list



#Verifica che il campo degli attributi sia corretto
#Ogni attributo deve terminare con un punto e virgola (;) ed essere separato dal successivo da esattamente uno spazio
#Ogni coppia nome_attributo e valore_attributo deve essere separata esattamente da uno spazio
#I valori non numerici devono essere inclusi in doppi apici
#attributi gene_id e transcript_id sono obbligatori e devono essere i primi due attributi di ogni riga
#Il valore di transcript_id nei record 'intron_CNS' deve essere diverso da ""
#Il valore di transcript_id nei record 'inter' e 'inter_CNS' deve essere uguale a ""
#Si sottolina anche che il nome degli attributi non dovrebbe essere incluso in doppi apici, infatti, 
#se ciò accadesse non verrebbero riconosciuti gli attributi obbligatori 
def check_attributes(row, errors):
    field_list = row[1].split('\t')
    mandatory_attributes = {'transcript_id': False, 'gene_id': False}   
    attribute_order = []
    attribute_name_with_dublequotes = False
    feature = field_list[2]

    attribute_list = re.findall(r'((?:(?:\"[^\"]*\")|\w)+(?:(?:\"[^\"]*\")|[\w\.]|\s)*)', field_list[8]) #Cerca ogni possibile coppia nome_attributo valore_attributo (esclude i separatori)
    attribute_list = [attribute.strip() for attribute in attribute_list]
    
    if '; '.join(attribute_list) + ';' != field_list[8]:   #controlla che alla fine di ogni attributo ci sia un punto e virgola e siano separati da uno spazio
        update_dict(errors[0], row[0], 'Illegal attribute separator: attributes must end in a semicolon which must then be separated by exactly one space character')
    else:
        for attribute in attribute_list:
            attribute_component = re.findall(r'((?:\"[^\"]*\")|[\w\.]+)', attribute)   #Divide ogni possibile coppia di attributi nei vari campi
            
            if len(attribute_component) != 2 or ' '.join(attribute_component) != attribute: #se ci sono + di due elementi o se non sono separati da esattamente uno spazio, dà errore 
                update_dict(errors[0], row[0], 'Attribute ' + attribute_component[0] + ' has the wrong format: each attribute must be a pair name value separated by exactly one space')
            elif not attribute_component[1].isnumeric() and (attribute_component[1][0] != '\"' or attribute_component[1][-1] != '\"'): #controlla che i valori non numerici siano inclusi in doppi apici
                update_dict(errors[0], row[0], 'Textual value of ' + attribute_component[0] + ' must be surrounded by doublequotes')
            
            if len(attribute_component) >= 2:  
                if attribute_component[0][0] == '\"' and attribute_component[0][-1] == '\"': #controlla se il nome dell'attributo è incluso in doppi apici
                    attribute_name_with_dublequotes = True

                if attribute_component[0] in mandatory_attributes:  #controlla se è un attributo obbligatorio
                    mandatory_attributes[attribute_component[0]] = True
                attribute_order.append(attribute_component[0])
                
                if attribute_component[0] == 'transcript_id':      
                    if feature == 'intron_CNS' and attribute_component[1] == '""':
                        update_dict(errors[0], row[0], "Illegal value of transcript_id: a 'intron_CNS' record musn't have transcript_id equal to \"\"")
                    if feature in ['inter_CNS', 'inter'] and attribute_component[1] != '""':
                        update_dict(errors[0], row[0], "Illegal value of transcript_id: a '" + feature + "' record must have transcript_id equal to \"\"")
        

        if len(attribute_order) >= 2 and (attribute_order[0] not in mandatory_attributes or attribute_order[1] not in mandatory_attributes): #controlla che i primi due attributi siano gene_id e transfert_id
            update_dict(errors[0], row[0], "transfert_id and gene_id must be the firsts two attributes")
        
        for key in mandatory_attributes:    #controlla che ci siano gli attributi obbligatori
            if not mandatory_attributes[key]:
                update_dict(errors[0], row[0],'Attribute ' + key + ' is required')
        if attribute_name_with_dublequotes:     #verifica se nella riga c'è almeno un nome di attributo incluso in doppi apici
            update_dict(errors[0], row[0], "WARNING: attribute names shouldn't be sourrounded by doublequotes")
    return                 



#Controlla che il valore del campo frame sia correto
#nel caso di start_codon e stop_codon non contigui il frame deve essere '0' o '1' o '2'
#nel caso siano contigui allora deve essere per forza 0
#nel caso deglle altre feature il valore del frame deve essere '0' o '1' o '2' o '.'
def check_frame(row, errors):
    frame = row[1].split('\t')[7]
    feature = row[1].split('\t')[2]
    if feature == 'start_codon' or feature == 'stop_codon':
        if row[0] not in errors[1]: #controlla che i valori di start e end della riga siano corretti per poter effettuare la conversione in intero
            start = int(row[1].split('\t')[3])
            end = int(row[1].split('\t')[4])
            if end - start > 0 and frame != '0':
                update_dict(errors[0], row[0], "Illegal frame (" + frame + "): contiguous " + feature + " must have frame '0'")
            elif end - start == 0 and frame not in ['0', '1', '2']:
                update_dict(errors[0], row[0], "Illegal frame (" + frame + "): " + feature + " must have a '0' or '1' or '2' in frame field")
    else:
        if frame not in ['0', '1', '2', '.']:
            update_dict(errors[0], row[0], "Illegal frame (" + frame + "): must have '0' or '1' or '2' or '3' or '.' in frame field")



#Verifica che il valore del campo strand sia corretto
#Deve essere per forza o '+' o '-'
def check_strand(row, errors):
    strand = row[1].split('\t')[6]
    if strand != '+' and strand != '-':
        update_dict(errors[0], row[0], "Illegal strand (" + strand + "): must have '+' or '-' in strand field")

#Verifica che il valore del campo score sia correto
#Deve essere o '.' o un intero o un floating point
def check_score(row, errors):
    score = row[1].split('\t')[5]
    if score != '.' and not re.match(r"([0-9]+(\.[0-9]+)?)|(\.[0-9]+)", score):
        update_dict(errors[0], row[0], "Illegal score (" + score + "): must have a number or a dot in score field")



#Verifiva la correttezza dei campi start e end
#Entrambi devono essere interi e maggiori o uguali di 1
#Inoltre deve risultare start <= end
def check_start_end(row, errors):
    start = row[1].split('\t')[3]
    end = row[1].split('\t')[4]
    line = row[0] 
    wrong_start_end = False
    if re.search("[^0-9]+", start):
        update_dict(errors[0], line, "Illegal start (" + start + "): must have an integer >= 1 in start field")
        wrong_start_end = True
    elif int(start) <= 0:
        update_dict(errors[0], line, "Illegal start (" + start + "): must have an integer >= 1 in start field")
        wrong_start_end = True

    if re.search("[^0-9]+", end):
        update_dict(errors[0], line, "Illegal end (" + end + "): must have an integer >= 1 in end field")
        wrong_start_end = True
    elif int(end) <= 0:
        update_dict(errors[0], line, "Illegal end (" + end + "): must have an integer >= 1 in end field")
        wrong_start_end = True

    if not wrong_start_end and int(start) > int(end):
        update_dict(errors[0], line, "Illegal value of start and end (" + start + ", " + end + "): it must be start <= end")
        wrong_start_end = True

    if wrong_start_end:
        errors[1].append(line)



#Controlla che i campi di ogni riga siano separati dal corretto separatore e che il numero di campi sia = 9
#inoltre controlla che siano presenti i record obbligatori 'CDS', 'start_codon' e 'stop_codon'
#verifica anche che la source nel file sia unica
#se questi primi vincoli sono rispettati, allora richiama per ogni riga le funzioni che verificano la correttezza dei vari campi che compongono la riga
#infine, ritorna un dizionario che ha come chiave la feature e
#come valore la lista di tutte le righe che hanno quella feature
#questo dizionario servirà per successivi controlli su vincoli legati a più record
def check_fields(rows, errors):
    required_feature = ['CDS', 'start_codon', 'stop_codon']
    row_dict = {'CDS': [], 'start_codon': [], 'stop_codon': [], '5UTR': [], '3UTR': [], 'inter': [], 'inter_CNS': [], 'intron_CNS': [], 'exon': []} #dizionario che verrà ritornato e che viene utilizzato per il controllo dei record obbligatori
    source_set = set()
    for row in rows:
        if re.search('\t', row[1]):         #verifica che sia presente il separatore dei campi all'interno della riga
            field_list = row[1].split('\t')     
            field_number = len(field_list)
            feature = field_list[2]
            if field_number != 9:           #verifica che il numero di campi sia esattamente 9
                update_dict(errors[0], row[0], "Wrong number of fields (" + str(field_number) + ") -> expected 9")
            elif feature in row_dict:       #se tutto ok e la feature è valida, richiama le funzioni per controllare i vari campi della riga
                    source_set.add(field_list[1])     
                    check_start_end(row, errors)
                    check_score(row, errors)
                    check_strand(row, errors)
                    check_frame(row, errors)
                    check_attributes(row, errors)
                    update_dict(row_dict, feature, row)
        else:
            update_dict(errors[0], row[0], "Illegal field separator")
    for feature in required_feature:
        if len(row_dict[feature]) == 0:
            update_dict(errors[0], 0, "The file needs a "+ feature +" record")
    if len(source_set) > 1:
        update_dict(errors[0], 0, "The source must be unique in the file")
    return row_dict



#verifica la somma delle lunghezze di tutti i record 'start_codon' siano <= 3bp
#controlla che le coordinate di tutti i record 'start_codon' siano all'interno delle coordinate di almeno un record 'CDS'
#inltre verifica che le coordinate di tutti i record 'start_codon' siano esterni alle coordinate di tutti i record '5UTR'
def check_start_codon(row_dict, errors):
    start_codon_length = 0
    start_codon_rows = []
    rows_not_in_cds = []    #contiene i numeri delle righe 'start_codon' che non sono incluse in coordinate di almeno un record 'CDS'
    rows_in_5UTR = []       #contiene i numeri delle righe che sono incluse nelle coordinate di almeno un record '5UTR'
    for row in row_dict['start_codon']:
        row_in_cds = False
        start_codon_rows.append(str(row[0]))
        if row[0] not in errors[1]: #controlla che i valori di start e end della riga siano corretti per poter effettuare la conversione in intero
            start = int(row[1].split('\t')[3])
            end = int(row[1].split('\t')[4])
            start_codon_length += end - start + 1
            for cds_row in row_dict['CDS']:
                if cds_row[0] not in errors[1]: #controlla che i valori di start e end della riga siano corretti per poter effettuare la conversione in intero
                    cds_start = int(cds_row[1].split('\t')[3])
                    cds_end = int(cds_row[1].split('\t')[4])
                    if cds_start <= start and end <= cds_end:
                        row_in_cds = True
                        break
            for row_5UTR in row_dict['5UTR']:   
                if row_5UTR[0] not in errors[1]: #controlla che i valori di start e end della riga siano corretti per poter effettuare la conversione in intero
                    start_5UTR = int(row_5UTR[1].split('\t')[3])
                    end_5UTR = int(row_5UTR[1].split('\t')[4])
                    if start in range(start_5UTR, end_5UTR + 1) or end in range(start_5UTR, end_5UTR + 1):
                        rows_in_5UTR.append(str(row[0]))     
                        break
            if not row_in_cds:  
                rows_not_in_cds.append(str(row[0]))
    if start_codon_length > 3:
        update_dict(errors[0], 0, 'The start_codon feature is up to 3bp long in total, violation at lines: ' + ', '.join(start_codon_rows))
    if len(rows_not_in_cds) > 0:
        update_dict(errors[0], 0, 'The start_codon must be included in coordinates for CDS features, violation at lines: ' + ', '.join(rows_not_in_cds))
    if len(rows_in_5UTR):
        update_dict(errors[0], 0, 'The stop codon must be exluded from the coordinates for the "5UTR" features, violation at lines: ' + ', '.join(rows_in_5UTR))
    return



#verifica la somma delle lunghezze di tutti i record 'stop_codon' siano <= 3bp
#controlla che le coordinate di tutti i record 'stop_codon'siano escluse delle coordinate di tutti i record 'CDS' e '3UTR'
def check_stop_codon(row_dict, errors):
    stop_codon_length = 0
    stop_codon_rows = []
    rows_in_cds = []    #contiene i numeri delle righe 'stop_codon' che sono incluse in coordinate di almeno un record 'CDS'
    rows_in_3UTR = []   #contiene i numeri delle righe che sono incluse nelle coordinate di almeno un record '3UTR'
    for row in row_dict['stop_codon']:
        stop_codon_rows.append(str(row[0]))
        if row[0] not in errors[1]: #controlla che i valori di start e end della riga siano corretti per poter effettuare la conversione in intero
            start = int(row[1].split('\t')[3])
            end = int(row[1].split('\t')[4])
            stop_codon_length += end - start + 1
            for cds_row in row_dict['CDS']:
                if cds_row[0] not in errors[1]: #controlla che i valori di start e end della riga siano corretti per poter effettuare la conversione in intero
                    cds_start = int(cds_row[1].split('\t')[3])
                    cds_end = int(cds_row[1].split('\t')[4])
                    if start in range(cds_start, cds_end + 1) or end in range(cds_start, cds_end + 1):
                        rows_in_cds.append(str(row[0]))
                        break
            for row_3UTR in row_dict['3UTR']:
                if row_3UTR[0] not in errors[1]: #controlla che i valori di start e end della riga siano corretti per poter effettuare la conversione in intero
                    start_3UTR = int(row_3UTR[1].split('\t')[3])
                    end_3UTR = int(row_3UTR[1].split('\t')[4])
                    if start in range(start_3UTR, end_3UTR + 1) or end in range(start_3UTR, end_3UTR + 1):
                        rows_in_3UTR.append(str(row[0]))
                        break
    if stop_codon_length > 3:
        update_dict(errors[0], 0, 'The stop_codon feature is up to 3bp long in total: violation at lines: ' + ', '.join(stop_codon_rows))
    if len(rows_in_cds) > 0:
        update_dict(errors[0], 0, 'The stop codon must not be included in the CDS features, violation at lines: ' + ', '.join(rows_in_cds))
    if len(rows_in_3UTR):
        update_dict(errors[0], 0, 'The stop codon must be exluded from the coordinates for the "3UTR" features, violation at lines: ' + ', '.join(rows_in_3UTR))
    return



#Si occupa di stampare le violazioni (se presenti)
def print_errors(file_input_name, errors):
    if errors == {}:
        print("The file '"+ file_input_name + "' is correct")
    else:
        print('Violations in file \'' + file_input_name + '\':')
        if 0 in errors:
            print("Gereral errors:")
            for error in errors[0]:
                print("\t" + error)
            del errors[0]
            print()
        for line in dict(sorted(errors.items())):
            print("At line " + str(line) + ":")
            for error in errors[line]:
                print("\t" + error)
            print()




#variabile composta da un dizionario e una lista
#il dizionario conterrà tutte le violazioni presenti: la chiave è il numero di riga e il valore è una lista che contiene tutte le violazioni di quella riga
#la lista conterrà tutte le righe che hanno valori irregolari di start o end o start > end, questa lista serve per verificare in modo rapido 
#se è possibile convertire in intero i campi start e end, nel caso questi valori servissero per qualche ulteriore controllo
errors = ({}, [])
file_input_name = input("File to validate: ")
file_rows = read_file(file_input_name)
format_file_rows = format_rows(file_rows)
row_dict = check_fields(format_file_rows, errors)
check_start_codon(row_dict, errors)
check_stop_codon(row_dict, errors)
print_errors(file_input_name, errors[0])