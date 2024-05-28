import sys
import json
import re

# This script convert the ngs_filter_file.json made by metabase api to the
# original ngsfilter file.
# arguments:
# -1: the file name of json file containing all information about tubes (libraries)
# -2: the output file name wihtout extension
# -3: the file name of tsv file containing pcr plate information


# parse the json file to get the reverse primer (get all reverse primer and
# get sequence if only one primer is found else return None)
def get_primer_reverse(json_file):
    with open(json_file) as jsf:
        data = json.load(jsf)
        primer_reverse_list = list()
        for plate in data['plates']:
            for well in plate["wells"]:
                if well['primer_reverse'] == 'undef':
                    continue
                if well['primer_reverse'] not in primer_reverse_list:
                    primer_reverse_list.append(well['primer_reverse'])
        if len(primer_reverse_list) > 1:
            return 'undef'
    return primer_reverse_list[0]


# main function: parse the json file and format data to the original ngsfilter file
def read_json_file(json_file, output_file):
    ofile = open(output_file + ".ngsfilter", "w")
    primer_reverse = get_primer_reverse(json_file)
    with open(json_file) as jsf:
        data = json.load(jsf)
        plate_number = 1
        rep_number = dict()
        plate_count = 0
        for plate in data['plates']:
            plate_count += 1
            plate_nb = str(plate_count).zfill(2)
            for well in plate["wells"]:
                if well["well_name"] == "extneg":
                    well['well_name'] = "coExt"
                if well['primer_reverse'] == 'undef':
                    well['primer_reverse'] = primer_reverse
                tag_combo = well['tag_forward'] + ":" + well['tag_reverse']
                position = ('F @ position=' + str(plate_nb) +
                            '_' + well["position"]["column"] +
                             well["position"]["row"] + ";")
                if well['well_type'] != 'sample':
                    if well['well_type'] == 'blank':
                        position += "type=control;control_type=sequencing;"
                    elif well['well_type'] == 'PCRPos':
                        position += "type=control;control_type=positive;"
                    elif well['well_type'] == 'PCRNeg':
                        position += "type=control;control_type=pcr;"
                    elif well['well_type'] == 'extneg':
                        position += "type=control;control_type=extraction;"
                else:
                    position += "type=sample;"
                if well['well_name'] not in rep_number:
                    rep_number[well['well_name']] = 1
                if well['well_type'] != 'sample':
                    extract = str(well['well_name']).replace(" ", "").replace("–", "_").replace('-', '_') + "_R" + str(rep_number[well['well_name']])
                else:
                    extract = str(well['well_name']).replace(" ", "").replace("–", "_").replace('-', '_')
                position += "sample_id=" + str(well['well_name']).replace(" ", "").replace("–", "_").replace('-', '_') + ";"
                rep_number[well['well_name']] += 1
                line = [
                    str(plate['name']).replace('-', '_'),
                    extract,
                    tag_combo,
                    well['primer_forward'],
                    well['primer_reverse'],
                    position
                ]
                ofile.write('\t'.join(line) + '\n')
            plate_number += 1
    ofile.close()

if __name__ == "__main__":
    read_json_file(sys.argv[1], sys.argv[2])
