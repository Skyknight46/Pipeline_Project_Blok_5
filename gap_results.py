import easygui
import taxonomy_script as tax
import re
import json


def fasta_lezen(fasta):
    """
    Fasta reader zet fasta in een dictionary met accessiecode als key
    en de sequentie als value. Als er hypothetical staat in de header
    wordt het niet meegenomen in de dictionary.
    :param fasta: fasta bestand
    :return: Dictionary met als key de volledige header en als value de sequentie
	en een lijst met alle accessiecodes die geen HYPOTHETICAL in de header hadden staan.
    """
    kijk = {}
    acc_code = []
    with open(fasta) as huidig:
        for line in huidig:
            if line.startswith(">") and "HYPOTHETICAL" not in line.upper():
                acc_code.append(line.split("/")[0].split(">")[1])
                key = line.strip()
                # print(key)
                kijk[key] = ""
            else:
                kijk[key] += (line.strip())

    print("*-*" * 5)
    print("Headers:", len(kijk.keys()))
    print("Sequenties:", len(kijk.values()))
    print("*-*" * 5)
    return kijk, acc_code


def gapcheck(kijk, acc_codes):
    """
    Kijkt in een aangegeven range (gap) of er aminozuren bevinden
    Deze
    :param kijk: Dictionary met als key de volledige header en als value de sequentie
    :param acc_codes: lijst met alle accessiecodes die geen HYPOTHETICAL in de header hadden staan
    :return: lijst met accessiecodes van sequenties met aminozuren in de gap
    """
    count = 0
    finalaccesies = []
    try:
        for key, value in kijk.items():
            position = value[901:916]
            # in iter 3 positie 901 - 915
            if not position == "-" * (916 - 901):
                if not acc_codes[count] in finalaccesies:
                    finalaccesies.append(acc_codes[count])
            count += 1
        print(len(finalaccesies))
        print("*-*" * 5)
    except IndexError:
        pass
    return finalaccesies


def taxonomy(finalaccesies):
    """
    Haalt de taxonomy op van finalaccessies en kijkt hoevaak
    eek bepaalde taxonomy voorkomt.
    :param finalaccesies: lijst met accessiecodes van sequenties met aminozuren in de gap
    """
    file_aantal = 0 # geeft aan bij welke file het script is.
    taxo = []
    for i in finalaccesies:
        data = tax.entrez_id_search(str(i))
        file_aantal += 1
        print("File ", file_aantal, " completed")

        for key, values in data.items():
            tax_list = []
            tax_re = re.search("(?<=\[).+?(?=\])", values)  # taxonomy
            tax_list.append(tax_re.group())
            taxo.append(tax_list)

    rijk = {}   # indeling taxonomy
    afdeling = {}
    klasse = {}
    orde = {}
    familie = {}
    geslacht = {}

    for j in taxo:
        try:
            if not j[0].split(",")[0].replace("'", "").strip(" ") in \
                   rijk.keys():
                rijk[j[0].split(",")[0].replace("'", "").strip(" ")] = 1
            else:
                rijk[j[0].split(",")[0].replace("'", "").strip(" ")] += 1
            if not j[0].split(",")[1].replace("'", "").strip(" ") in \
                   afdeling.keys():
                afdeling[j[0].split(",")[1].replace("'", "").strip(" ")] = 1
            else:
                afdeling[j[0].split(",")[1].replace("'", "").strip(" ")] += 1
            if not j[0].split(",")[2].replace("'", "").strip(" ") in \
                   klasse.keys():
                klasse[j[0].split(",")[2].replace("'", "").strip(" ")] = 1
            else:
                klasse[j[0].split(",")[2].replace("'", "").strip(" ")] += 1
            if not j[0].split(",")[3].replace("'", "").strip(" ") in \
                   orde.keys():
                orde[j[0].split(",")[3].replace("'", "").strip(" ")] = 1
            else:
                orde[j[0].split(",")[3].replace("'", "").strip(" ")] += 1
            if not j[0].split(",")[4].replace("'", "").strip(" ") in \
                   familie.keys():
                familie[j[0].split(",")[4].replace("'", "").strip(" ")] = 1
            else:
                familie[j[0].split(",")[4].replace("'", "").strip(" ")] += 1
            if not j[0].split(",")[5].replace("'", "").strip(" ") in \
                   geslacht.keys():
                geslacht[j[0].split(",")[5].replace("'", "").strip(" ")] = 1
            else:
                geslacht[j[0].split(",")[5].replace("'", "").strip(" ")] += 1
        except IndexError:
            pass
    print(rijk, " Rijk")
    print(afdeling, " Afdeling")
    print(klasse, " Klasse")
    print(orde, " Orde")
    print(familie, " Familie")
    print(geslacht, " Geslacht")
    bestand = "AnalyseTaxonomy.txt"
    with open(bestand, 'w') as file:
        file.write(json.dumps(rijk) + " Rijk" + "\n")
        file.write(json.dumps(afdeling) + " Afdeling" + "\n")
        file.write(json.dumps(klasse) + " Klasse" + "\n")
        file.write(json.dumps(orde) + " Orde" + "\n")
        file.write(json.dumps(familie) + " Familie" + "\n")
        file.write(json.dumps(geslacht) + " Geslacht" + "\n")


def main():
    fasta = easygui.fileopenbox()
    kijk, volledige_header = fasta_lezen(fasta)
    finalaccesies = gapcheck(kijk, volledige_header)
    taxonomy(finalaccesies)


main()
