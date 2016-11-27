import csv

# merge 2 .csv files into single

semiglo = []
local = []

def readCSV(filename, destination) :
  print(filename)
  with open(filename, 'r') as csvfile : 
    filereader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in filereader :
      destination.append(row)
  print(destination[0])

def write(first, second) :
  with open('merged.csv', 'w') as csvfile:
      filewriter = csv.writer(csvfile, delimiter=' ',
                              quotechar='|', quoting=csv.QUOTE_MINIMAL)

      for i in range(len(first)) :
        filewriter.writerow(first[i])
        filewriter.writerow(second[i]);
        filewriter.writerow([]);


lfile = 'Cas_Local_Alignment.csv'
smfile = 'Cas_SemiGlobal_Alignment.csv'


readCSV(smfile, semiglo)
readCSV(lfile, local)
write(semiglo, local)

