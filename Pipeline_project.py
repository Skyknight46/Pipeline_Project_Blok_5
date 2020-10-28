# Pipeline project course 5, project groep 6
# Zorg dat het script gedraaid word vanuit de map waar ook het script esl-reformat 
# in te vinden is. standaard is dit: ....\dir\easel\miniapps\esl-reformat.c
import os
import sys
import subprocess
import FastaEditor as fasta_add
from sys import argv


def clustalo(clustalo_input_file, clustalo_output_file):
    """
    Zorgt voor het command om ClustalO te draaien.
    :param clustalo_input_file: pad naar input file
    :param clustalo_output_file: pad naar output file
    :return: pad naar MSA van inputsequenties
    """
    cmd = "clustalo -o {} -i {} -v".format(clustalo_output_file,
                                           clustalo_input_file)
    e = subprocess.check_call(cmd, shell=True)
    return clustalo_output_file


def hmm_build(clustalo_output_file, hmmbuild_output, cpucores):
    """
    Maakt met de output van ClustalO een hmmbuild aan.
    :param clustalo_output_file: pad naar MSA met inputsequenties
    :param hmmbuild_output: pad naar output file (HMM logo)
    :param cpucores: aantal threads van pc -2
    :return: filepath van de outputfile
    """
    cmd = "hmmbuild --amino --cpu {} {} {}".format(cpucores, hmmbuild_output,
                                                   clustalo_output_file)
    e = subprocess.check_call(cmd, shell=True)
    return hmmbuild_output


def hmm_search(counter, stockholm_file_hmmsearch, hmmbuild_output, refseq_pwd,
               esl_format_fasta, cpucores):
    """
    Hmmsearch zoekt met het HMM-logo van HMMbuild tegen de RefSeq_db naar homologe 
	sequenties. 
	esl-reformat zorgt ervoor dat de stockholm output van hmmsearch wordt omgezet 
	naar fasta.
    :param counter: counter om de unieke naamgeving van files te verzorgen.
    :param stockholm_file_hmmsearch: pad naar msafile met stockholm bestandsindeling
    :param hmmbuild_output: pad naar HMM-logo file (Output van HMMbuild)
    :param refseq_pwd: pad naar de Refseq db (single fasta file)
    :param esl_format_fasta: pad van de naar fasta geconverteerde stockholm msa file
    :param cpucores: aantal threads van pc -2
	:return: pad naar file met outputsequenties in fasta format
    """

    print("-=" * 30)
    print("hmmsearch")
    print("stockholm_file_hmmsearch: ", stockholm_file_hmmsearch)
    print("hmmbuild_output: ", hmmbuild_output)
    print("refseq_pwd:", refseq_pwd)
    print("aantal cores: " + str(cpucores))
    print("-=" * 30)
    outputbestand = "~/main_output_hmmsearch" + str(counter) + ".txt"
    cmd_hmmsearch = "hmmsearch -o {} -A {} --notextw --acc --noali " \
                    "--cpu {} {} {}".format(outputbestand,
                                            stockholm_file_hmmsearch, cpucores,
                                            hmmbuild_output, refseq_pwd)
    f = subprocess.check_call(cmd_hmmsearch, shell=True)
    cmd_remove_testbestand = "rm {}".format(outputbestand)
    g = subprocess.check_call(cmd_remove_testbestand, shell=True)
    print("#esl-reformat")
    cmd_format = "./esl-reformat -o {} --informat stockholm fasta " \
                 "{}".format(esl_format_fasta, stockholm_file_hmmsearch)
    e = subprocess.check_call(cmd_format, shell=True)
    return esl_format_fasta


def fasta_joiner(clustalo_input, esl_reformat, aanpasbare_file_location):
    """
    Roept de functie fasta_adder aan en voegt het orginele bestand en
    de output van esl-reformat samen.
    :param clustalo_input: pad naar de input_file van de huidige iteratie
    :param esl_reformat: pad naar de file met output van de huidige iteratie
    :param aanpasbare_file_location: pad naar de file die bij elke iteratie 
	overschreven wordt
    :return: pad naar de file die bij elke iteratie	overschreven wordt
    """
    if os.path.isfile(
            aanpasbare_file_location + "/aanpasbare_sequentie_file.fasta"):
        aanpasbare_file_location = aanpasbare_file_location + \
                                   "/aanpasbare_sequentie_file.fasta"
        overschrijfbare_file = fasta_add.fasta_adder(clustalo_input,
                                                     esl_reformat,
                                                     aanpasbare_file_location)

    elif not os.path.isfile(
            aanpasbare_file_location + "/aanpasbare_sequentie_file.fasta"):
        cmd_file_aanmaken = "touch {}/aanpasbare_sequentie_file.fasta".format(
            aanpasbare_file_location)
        subprocess.check_call(cmd_file_aanmaken, shell=True)
        aanpasbare_file_location = aanpasbare_file_location + \
                                   "/aanpasbare_sequentie_file.fasta"
        overschrijfbare_file = fasta_add.fasta_adder(clustalo_input,
                                                     esl_reformat,
                                                     aanpasbare_file_location)
    return overschrijfbare_file


def main():
    print("\n" * 2)
    pwd_check = os.getcwd()
    cpucores = os.cpu_count() - 2
    aanpasbare_file_location = argv[7]
    # hier alleen de path naar de opslaglocatie naar keuze zetten,
    # niet ook nog de bestandsnaam!!  ps. mocht je pc bij
    # het automatisch aanvullen van de map er een "/"
    # achteraanplakken... deze dan even weghalen.
    # Voorbeeld: <~/Desktop> of <~/Downloads/project_blok_5> dus niet
    # <~/Desktop/aanpasbarefile.txt> of <~/Desktop/>
    if "easel/miniapps" in pwd_check:
        pass
    else:
        print("#Je zit niet in de juiste map!\n#Huidige map: "
              + pwd_check + "\n#Als je easel nog niet hebt geïnstaleerd "
                            "download het via de volgende link:\n\t\t
							 https://github.com/EddyRivasLab/easel")
        sys.exit(-1)
    counter = 0
    if counter == 0:
        try:
            # eerst de naam van je output vervolgens naam bestand
            clustalo_input = argv[2] + ".txt"
            clustalo_output = argv[1] + str(counter) + ".txt"
            # print("#" + clustalo_output)
            hmmbuild_output = argv[3] + str(counter) + ".txt"  # naam
            stockholm_file_hmmsearch = argv[4] + str(counter) + ".txt"
            refseq_pwd = argv[5]  # pwd naar refseq
            esl_format_fasta = argv[6] + str(counter) + ".fasta"  # naam

            if not os.path.isfile(clustalo_input):
                print("#Input file ", clustalo_input, " not found!")
                print("\n")
                sys.exit(-1)
                # controleren of clustal input er al is
            else:
                pass

            if not os.path.isfile(clustalo_output):
                print("#" + clustalo_output, " not found \nrunning Clustalo.")
                print("-=-" * 20)
                clustalo(clustalo_input, clustalo_output)
                print("-=-" * 20)
                print("\n")

            else:
                print("#Found file: ", clustalo_output, "\n")

            if not os.path.isfile(hmmbuild_output):
                hmm_build(clustalo_output, hmmbuild_output, cpucores)
                print("\n")

            else:
                pass

            if not os.path.isfile(esl_format_fasta):
                esl_reformat = hmm_search(counter, stockholm_file_hmmsearch,
                                          hmmbuild_output, refseq_pwd,
                                          esl_format_fasta, cpucores)
                print("\n")

            else:
                pass
        except IndexError:
            print("Incorrect number of arguments specified\n"
                  "Please check if the arguments were correct")
            pass

        except:  # Vangt onbekende exceptions af
            print("#Unkown error, please try again")
            sys.exit(-1)

    clustalo_input = fasta_joiner(clustalo_input, esl_reformat,
                                  aanpasbare_file_location)
    counter += 1

    while counter != 3:  # itereert over de functies heen
        clustalo_output = argv[1] + str(counter) + ".txt"
        print("#" + clustalo_output)
        hmmbuild_output = argv[3] + str(counter) + ".txt"
        stockholm_file_hmmsearch = argv[4] + str(counter) + ".txt"
        esl_format_fasta = argv[6] + str(counter) + ".fasta"

        if not os.path.isfile(clustalo_output):
            print("#" + clustalo_output, " not found \n#running Clustalo")
            # PLACEHOLDER clustalo_input = pad naar clustalo_inputfile
            clustalo(clustalo_input, clustalo_output)
            print("\n")

        else:
            print("#Found file: ", clustalo_output, "\n")

        if not os.path.isfile(hmmbuild_output):
            hmm_build(clustalo_output, hmmbuild_output, cpucores)
            print("\n")

        else:
            pass

        if not os.path.isfile(esl_format_fasta):
            esl_reformat = hmm_search(counter, stockholm_file_hmmsearch,
                                      hmmbuild_output, refseq_pwd,
                                      esl_format_fasta, cpucores)
            print("\n")
            clustalo_input = fasta_joiner(clustalo_input, esl_reformat,
                                          aanpasbare_file_location)

        else:
            pass

        counter += 1

    #  na laatste iter runnen clustalo.
    clustalo(clustalo_input, "~/Desktop/clustal_out2_10.txt")


main()
