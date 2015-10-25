#! /usr/bin/python3
#-*-coding:utf-8-*-
#####DNA sequences region finding with sam file.


class Genome:  #generate the genome information
    def __init__(self, filename):
        self.file = open(filename)
        self.header = self.file.readline().rstrip().split(sep = ',')  #header list
        print("The file's header is", self.header, sep = ' ')
        self.__headerlen = len(self.header)
        self.chr = {}
        for line in self.file:
            self.chr[line.rstrip().split(sep = ',')[0]] = line.rstrip().split(sep = ',')[1:-1]
    def chromosomes(self): #return the chromosome names
        return self.chr.keys()
    def length(self, chr, length = 'Length(mm)'):  #find the length of chromosome
        return float(self.chr[chr][0])
    def basepairs(self, chr, basepairs = 'Basepairs'):  #find the basepairs of chromosome
        return int(self.chr[chr][1])

class Distribution:  #bin distribution data
    def __init__(self, file, bin = 1000):
        self.__bin = bin
        self.genome = Genome(file)
        self.chrname = self.genome.chromosomes()
        #print(self.chrname)
        self.table = {}
        for chrname in self.chrname:  #initiate the table
            num = int(int(self.genome.basepairs(chrname)) / self.__bin)
            #print(num)
            self.table[chrname] = []
            for i in range(0, num+1):
                self.table[chrname].append(0)  #i should be converted to char
    def calculate(self, inputfile):  #calculate the result
        self.calinput = open(inputfile)
        #self.tempoutput = open('tempout.csv', 'w+')
        for line in self.calinput:
            if line[0] == '@': continue
            linechrname = line.split(sep = '\t')[2]  #chrname from samfile
            #print(linechrname, type(linechrname), end = '\r')
            if linechrname not in self.chrname: continue
            linestart = line.split(sep = '\t')[3]
            #print(linestart, end = '\r')
            linenum = int(int(linestart)/self.__bin)
            self.table[linechrname][linenum] += 1
            #print(self.table[linechrname][linenum])
        self.calinput.close()
        print('Calculate over!')
    def writefile(self, outputfile):  #write the resul
        self.caloutput = open(outputfile, 'w+')
        for chrname in self.chrname:
            num = int(int(self.genome.basepairs(chrname)) / self.__bin)
            for i in range(0, num+1):
                if self.table[chrname][i] == 0: continue
                self.caloutput.write(chrname)
                self.caloutput.write(',')
                self.caloutput.write(str(i * self.__bin))
                self.caloutput.write(',')
                self.caloutput.write(str((i+1) * self.__bin))
                self.caloutput.write(',')
                self.caloutput.write(str(self.table[chrname][i]))
                self.caloutput.write('\n')
        self.caloutput.close()
        print('Result has been writed to', outputfile, sep = ' ')
    

#EOF
 
        
        
        
