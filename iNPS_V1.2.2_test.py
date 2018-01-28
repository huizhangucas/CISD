#!/usr/bin/env python
#--coding:utf-8 --


#################################################################################################
#################################################################################################
########                                                                                 ########
########    An improved Algorithm for Nucleosome Positioning from Sequencing Data        ########
########                                                                                 ########
########    Program to position nucleosomes with Keji Zhao's tag coodinate bed files.    ########
########                                                                                 ########
########    Author:  CHEN Weizhong                                                       ########
########                                                                                 ########
########    Working Environment:  Python3                                                ########
########                                                                                 ########
########    Date:      2012-05-02                                                        ########
########                                                                                 ########
########    Modified:  2013-11-15                                                        ########
########                                                                                 ########
########    Modified:  2014-02-24  Version 1.0 , Add Poisson test                        ########
########                                                                                 ########
########    Update:    2014-10-25  Version 1.1                                           ########
########                                                                                 ########
########    Update:    2014-11-13  Version 1.1.2                                         ########
########                                                                                 ########
########    Update:    2015-01-16  Version 1.2.0 , Available for paired-end data         ########
########                                                                                 ########
########    Update:    2015-03-03  Version 1.2.1                                         ########
########                                                                                 ########
########    Update:    2015-08-21  Version 1.2.2                                         ########
########                                                                                 ########
#################################################################################################
#################################################################################################








### Input and output files parameters
FilesParameters={'inputfile_for_nucleosome_positioning':'none',
                 'Chosen_Chromosome_Abbreviation':'none',
                 'Chosen_Chromosome_Length':'none',
                 'outputfile_like_bed':'',
                 'outputfile_like_wig':'',
                 }




### The standard deviation in Gaussian convolution ### Please set them with integers.
PreliminaryConvolutionParameters={'sigma':3,
                                  'times_of_sigma':3,
                                  'secondary_sigma':1}

### Threshold of height and width of peak
threshold={'PCC_independence_1':0.8,   'concavity_and_convexity_1':0.8,  'gap_length_1':100,  'height_ratio_1':0,
           'PCC_independence_2':0.2,   'concavity_and_convexity_2':0.8,  'gap_length_2':10,   'height_ratio_2':0.2,
           'PCC_independence_3':-0.5,  'concavity_and_convexity_3':0.8,  'gap_length_3':6,    'height_ratio_3':0.4,
           'PCC_independence_4':-1,    'concavity_and_convexity_4':0.8,  'gap_length_4':20,   'height_ratio_4':0.5,
           'PCC_independence_5':-1,    'concavity_and_convexity_5':0,    'gap_length_5':6,    'height_ratio_5':0.7,
           'PCC_independence_6':-1,    'concavity_and_convexity_6':0,    'gap_length_6':10,   'height_ratio_6':0.8,
           ########
           'influence_coefficient_cutoff':0,  'influence_coefficient_absolute_cutoff':0,
           'influence_coefficient_distance':9,  'influence_coefficient_distance_center':13,
           ########
           'filter_switch':'on',  'merging_center_distance':150,  'merging_height_watershed':5,
           'merging_height_ratio_1':0.80,  'merging_gap_1':50,  'merging_percentage_1':0.75,
           'merging_height_ratio_2':0.80,  'merging_gap_2':70,  'merging_percentage_2':0.75,
           'merging_height_ratio_3':0.90,  'merging_gap_3':70,  'merging_percentage_3':0.60,
           'merging_height_ratio_4':0.60,  'merging_gap_4':40,  'merging_percentage_4':0.70,
           'discarded_noise_selfAUC_1':11,    'discarded_noise_selfLength_1':5,   'discarded_noise_LoG_average_1':-0.06,
           'discarded_noise_selfLength_2':5,  'discarded_noise_LoG_average_2':-0.03,
           'discarded_noise_selfHeight_3':3,  'discarded_noise_LoG_average_3':-0.01,
           'discarded_noise_RealSelfLength_4':2,  'discarded_noise_LoG_average_4':-100.0,
           'discarded_noise_RealSelfLength_5':3,  'discarded_noise_LoG_average_5':-0.5,
           'discarded_noise_RealSelfLength_6':4,  'discarded_noise_LoG_average_6':-0.15}

########chromosome length in hg18
########chromosome_length={'chr1':247249719,'chr2':242951149,'chr3':199501827,'chr4':191273063,'chr5':180857866,'chr6':170899992,'chr7':158821424,'chr8':146274826,'chr9':140273252,'chr10':135374737,'chr11':134452384,'chr12':132349534,'chr13':114142980,'chr14':106368585,'chr15':100338915,'chr16':88827254,'chr17':78774742,'chr18':76117153,'chr19':63811651,'chr20':62435964,'chr21':46944323,'chr22':49691432,'chrX':154913754,'chrY':57772954}










### Main program
import math
import time
import sys
from optparse import OptionParser
import os

def Gaussian_profile(PreliminaryConvolutionParameters):
    sigma=PreliminaryConvolutionParameters['sigma']
    times_of_sigma=PreliminaryConvolutionParameters['times_of_sigma']
    secondary_sigma=PreliminaryConvolutionParameters['secondary_sigma']
    Gaussian=[]
    First_Derivative_of_Gaussian=[]
    LoG=[]
    Third_Derivative_of_Gaussian=[]
    for x0 in range(-sigma*times_of_sigma,sigma*times_of_sigma+1):
        x=x0*(-1)
        Gaussian.append(math.exp(-x*x/(2.0*sigma*sigma))/math.sqrt(2*(math.pi)*sigma*sigma))
        First_Derivative_of_Gaussian.append(((-x)/((sigma*sigma)*(math.sqrt(2*(math.pi)*sigma*sigma))))*(math.exp(-x*x/(2.0*sigma*sigma))))
        LoG.append(((x*x)/(sigma*sigma*sigma*sigma)-1/(sigma*sigma))*(math.exp(-x*x/(2.0*sigma*sigma)))/(math.sqrt(2*(math.pi)*sigma*sigma)))
        Third_Derivative_of_Gaussian.append(((3*x)/(sigma*sigma*sigma*sigma)-(x*x*x)/(sigma*sigma*sigma*sigma*sigma*sigma))*(math.exp(-x*x/(2.0*sigma*sigma)))/(math.sqrt(2*(math.pi)*sigma*sigma)))
    secondary_Gaussian=[]
    secondary_First_Derivative_of_Gaussian=[]
    secondary_LoG=[]
    secondary_Third_Derivative_of_Gaussian=[]
    for x0 in range(-secondary_sigma*times_of_sigma,secondary_sigma*times_of_sigma+1):
        x=x0*(-1)
        secondary_Gaussian.append(math.exp(-x*x/(2.0*secondary_sigma*secondary_sigma))/math.sqrt(2*(math.pi)*secondary_sigma*secondary_sigma))
        secondary_First_Derivative_of_Gaussian.append(((-x)/((secondary_sigma*secondary_sigma)*(math.sqrt(2*(math.pi)*secondary_sigma*secondary_sigma))))*(math.exp(-x*x/(2.0*secondary_sigma*secondary_sigma))))
        secondary_LoG.append(((x*x)/(secondary_sigma*secondary_sigma*secondary_sigma*secondary_sigma)-1/(secondary_sigma*secondary_sigma))*(math.exp(-x*x/(2.0*secondary_sigma*secondary_sigma)))/(math.sqrt(2*(math.pi)*secondary_sigma*secondary_sigma)))
        secondary_Third_Derivative_of_Gaussian.append(((3*x)/(secondary_sigma*secondary_sigma*secondary_sigma*secondary_sigma)-(x*x*x)/(secondary_sigma*secondary_sigma*secondary_sigma*secondary_sigma*secondary_sigma*secondary_sigma))*(math.exp(-x*x/(2.0*secondary_sigma*secondary_sigma)))/(math.sqrt(2*(math.pi)*secondary_sigma*secondary_sigma)))
    ConvolutionParameters={'Gaussian':Gaussian,
                           'First_Derivative_of_Gaussian':First_Derivative_of_Gaussian,
                           'LoG':LoG,
                           'Third_Derivative_of_Gaussian':Third_Derivative_of_Gaussian,
                           'sigma':sigma,
                           'times_of_sigma':times_of_sigma,
                           'secondary_Gaussian':secondary_Gaussian,
                           'secondary_First_Derivative_of_Gaussian':secondary_First_Derivative_of_Gaussian,
                           'secondary_LoG':secondary_LoG,
                           'secondary_Third_Derivative_of_Gaussian':secondary_Third_Derivative_of_Gaussian,
                           'secondary_sigma':secondary_sigma}
    return ConvolutionParameters





def log10( x ):
    return math.log(x)/math.log(10)

class Poisson_test:
    def greater( background , X ):
        if X > 0:
            Added = math.exp( - background )
            #Added = 1
            for i in range( 1 , X ):
                Added = Added * background /i
            cdf = 0
            i = X
            while Added > 0:
                Added = Added * background / i
                cdf = cdf + Added
                i = i + 1
            #cdf = cdf * math.exp( - background )
            score = math.log(cdf)/math.log(10)
            score = round( -score , 8 )
        elif X == 0:
            cdf = 1
            score = 0
        return cdf , score

    def less( background , X ):
        Added = math.exp( - background )
        #Added = 1
        cdf = 0 + Added
        for i in range( 1 , X+1 ):
            Added = Added * background / i
            cdf = cdf + Added
        #cdf = cdf * math.exp( -1 background )
        score = math.log(cdf)/math.log(10)
        score = round( -score , 8 )
        return cdf , score

    def greater_fast( background , X ):
        if X > 0:
            using = 100
            resting = background - using
            Added = math.exp( - using )
            for i in range( 1 , X ):
                Added = Added * background /i
                if Added > 0:
                    while log10( Added ) > 200:
                        Added = Added * math.exp(-10)
                        resting = resting - 10
                    if log10( Added ) < 10:
                        Added = Added * math.exp(10)
                        resting = resting + 10
                else:
                    pass
            cdf = 0
            i = X
            while Added > 0:
                Added = Added * background / i
                cdf = cdf + Added
                i = i + 1
                if Added > 0:
                    while log10( Added ) > 200:
                        Added = Added * math.exp(-100)
                        cdf = cdf * math.exp(-100)
                        resting = resting - 100
                else:
                    pass
            score = log10( cdf )
            cdf = cdf * math.exp( - resting )
            score = score - resting / math.log(10)
            score = round( -score , 8 )
        elif X == 0:
            cdf = 1
            score = 0
        return cdf , score

    def less_fast( background , X ):
        using = 700
        resting = background - using
        Added = math.exp( - using )
        cdf = 0 + Added
        for i in range( 1 , X+1 ):
            Added = Added * background / i
            cdf = cdf + Added
            if Added > 0:
                while log10( Added ) > 200:
                    Added = Added * math.exp(-100)
                    cdf = cdf * math.exp(-100)
                    resting = resting - 100
            else:
                pass
        score = log10( cdf )
        cdf = cdf * math.exp( - resting )
        score = score - resting / math.log(10)
        score = round( -score , 8 )
        return cdf , score







class MainProgram:
    def __init__(self,FilesParameters,ConvolutionParameters,threshold):
        print(time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        self.inputfile=FilesParameters['inputfile_for_nucleosome_positioning']
        self.celltype=self.inputfile.split('/')[-1].split('Nucleosome')[0]
        self.chromosome=self.inputfile.split('/')[-1].split('.')[0].split('-')[-1]
        self.chrlength=chromosome_length[self.chromosome]
        self.chrRecordingListLength=0
        self.outputfile_path=FilesParameters['fileout_path']
        self.outputfile_wiggle=self.outputfile_path+self.inputfile.split('/')[-1][0]+'_'+self.chromosome+FilesParameters['outputfile_suffix']+'.like_wig'
        self.outputfile_nucleosome=self.outputfile_path+self.inputfile.split('/')[-1][0]+'_'+self.chromosome+FilesParameters['outputfile_suffix']+'.like_bed'
        print('Input file:  ',self.inputfile)
        print('Cell type:   ',self.celltype+' CD4+ T cells from human')
        print('Chromosome:  ',self.chromosome)
        print('Length of',self.chromosome,'in hg18:  ',self.chrlength)
        print('Output wiggle file of nucleosome distribution: ',self.outputfile_wiggle)
        print('Output file of nucleosome collection:          ',self.outputfile_nucleosome)
        self.ConvolutionParameters=ConvolutionParameters
        self.threshold=threshold
        self.score_list=[]
        self.tag_list = [ ]
        self.Gaussian_list=[]
        self.FDoG_list=[]
        self.LoG_list=[]
        self.TDoG_list=[]
        self.score_table=[]    ### coordinate, self.score_list, self.Gaussian_list, self.LoG_list, column for recording
        ##self.secondary_Gaussian_list=[]
        ##self.secondary_FDoG_list=[]
        self.secondary_LoG_list=[]
        ##self.secondary_TDoG_list=[]
        print('======> ======> ======> ======> ======> ======>')


    def score(self):
        print('Reading tag coodinate bed file and scoring the nucleosomes ......')
        if self.chrlength%10==0:
            self.chrRecordingListLength=self.chrlength//10
            self.score_list=[0]*(self.chrlength//10)
            self.tag_list = [0]*( self.chrlength//10 )
        elif self.chrlength%10!=0:
            self.chrRecordingListLength=self.chrlength//10+1
            self.score_list=[0]*(self.chrlength//10+1)
            self.tag_list = [0]*( self.chrlength//10+1 )
        line=0
        text=open(self.inputfile,'r')
        while True:
            try:
                data=next(text).split('\n')[0].split('\r')[0].split('\t')
                line=line+1
#                sys.stdout.write('\r        Line:'+str(line) )

                if data[5]=='+':
                    #beginning=int(data[1])+37
                    #ending=int(data[1])+37+75
                    #TagCoordinate = int( data[1] ) + 75
                    TagCoordinate = int( (int(data[1])+int(data[2]))/2 ) + 58
                    beginning=TagCoordinate-37
                    ending=TagCoordinate+37
                elif data[5]=='-':
                    #beginning=int(data[2])-37-75
                    #ending=int(data[2])-37
                    #TagCoordinate = int( data[2] ) - 75
                    TagCoordinate = int( (int(data[1])+int(data[2]))/2 ) - 58
                    beginning=TagCoordinate-37
                    ending=TagCoordinate+37

                ####
                for i in range(beginning,ending+1):
                    if ( ((i-1)//10)<=self.chrRecordingListLength-1 ) and ( ((i-1)//10)>=0 ):
                        self.score_list[(i-1)//10]+=0.1
                    else:
                        ####print('Coordinate error in ',self.inputfile,'line:',line,' ',data, i, self.chrRecordingListLength-1)      #### #### #### #### 20141021 #### #### #### ####
                        pass                                                                                                          #### #### #### #### 20141021 #### #### #### ####
                ####
                if ( ((TagCoordinate-1)//10) <= self.chrRecordingListLength - 1 ) and ( ((TagCoordinate-1)//10)>=0 ):
                    self.tag_list[ (TagCoordinate-1)//10 ] += 1
                else:
                    ####print('Coordinate error in ',self.inputfile,'line:',line,' ',data, i, self.chrRecordingListLength-1)
                    pass
                ####
                sys.stdout.flush()
            except StopIteration:
                break
        text.close()
        print('...... File reading and nucleosome scoring is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        print('\tChecking ==> Total tag fragments:  ',line)
        print('\tChecking ==> Length of original self.score_list:  ',len(self.score_list))
        print('\tChecking ==> Length of original self.tag_list:    ',len(self.tag_list) )


    def score_for_Paired_end(self):
        print('Read paired-end tags and scoring the nucleosomes ......')
        if self.chrlength%10==0:
            self.chrRecordingListLength=self.chrlength//10
            self.score_list=[0]*(self.chrlength//10)
            self.tag_list = [0]*( self.chrlength//10 )
        elif self.chrlength%10!=0:
            self.chrRecordingListLength=self.chrlength//10+1
            self.score_list=[0]*(self.chrlength//10+1)
            self.tag_list = [0]*( self.chrlength//10+1 )
        line=0
        text=open(self.inputfile,'r')
        while True:
            try:
                data=next(text).split('\n')[0].split('\r')[0].split('\t')
                line=line+1
#                sys.stdout.write('\r        Line:'+str(line) )
                #if data[5]=='+':
                #    beginning=int(data[1])+37
                #    ending=int(data[1])+37+75
                #    TagCoordinate = int( data[1] ) + 75
                #elif data[5]=='-':
                #    beginning=int(data[2])-37-75
                #    ending=int(data[2])-37
                #    TagCoordinate = int( data[2] ) - 75
                A = int(data[1])
                B = int(data[2])
                tag_length = B - A + 1
                if ( tag_length >= FilesParameters[ 'pe_min' ] ) and ( tag_length <= FilesParameters[ 'pe_max' ] ):
                    beginning = round( A + 0.25*(B-A) )
                    ending    = round( B - 0.25*(B-A) )
                    TagCoordinate = round( (A+B)*0.5 )
                    ####
                    for i in range(beginning,ending+1):
                        if ( ((i-1)//10)<=self.chrRecordingListLength-1 ) and ( ((i-1)//10)>=0 ):
                            self.score_list[(i-1)//10]+=0.1
                        else:
                            ####print('Coordinate error in ',self.inputfile,'line:',line,' ',data, i, self.chrRecordingListLength-1)      #### #### #### #### 20141021 #### #### #### ####
                            pass                                                                                                          #### #### #### #### 20141021 #### #### #### ####
                    ####
                    if ( ((TagCoordinate-1)//10) <= self.chrRecordingListLength - 1 ) and ( ((TagCoordinate-1)//10)>=0 ):
                        self.tag_list[ (TagCoordinate-1)//10 ] += 1
                    else:
                        ####print('Coordinate error in ',self.inputfile,'line:',line,' ',data, i, self.chrRecordingListLength-1)
                        pass
                    ####
                else:
                    pass
                sys.stdout.flush()
            except StopIteration:
                break
        text.close()
        print('...... File reading and nucleosome scoring is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        print('\tChecking ==> Total tag fragments:  ',line)
        print('\tChecking ==> Length of original self.score_list:  ',len(self.score_list))
        print('\tChecking ==> Length of original self.tag_list:    ',len(self.tag_list) )
        


    def Gaussian_convolution_smoothing(self):
        print('Performing Gaussian convolution ......')
        def convolution(x,Gaussian,First_Derivative_of_Gaussian,LoG,Third_Derivative_of_Gaussian,sigma,times_of_sigma):    ### x is the index number of self.score_list, which is transmitted by parameter j in the following function "Gaussian_convolution_smoothing".
            y=0
            FDoG_y=0
            LoG_y=0
            TDoG_y=0
            if x>=sigma*times_of_sigma and x<=len(self.score_list)-sigma*times_of_sigma-1:
                for n in range(x-sigma*times_of_sigma,x+sigma*times_of_sigma+1):
                    y=y+self.score_list[n]*Gaussian[n-(x-sigma*times_of_sigma)]
                    FDoG_y=FDoG_y+self.score_list[n]*First_Derivative_of_Gaussian[n-(x-sigma*times_of_sigma)]
                    LoG_y=LoG_y+self.score_list[n]*LoG[n-(x-sigma*times_of_sigma)]
                    TDoG_y=TDoG_y+self.score_list[n]*Third_Derivative_of_Gaussian[n-(x-sigma*times_of_sigma)]
            elif x<sigma*times_of_sigma:
                for n in range(0,x+sigma*times_of_sigma+1):
                    y=y+self.score_list[n]*Gaussian[sigma*times_of_sigma-x+n]
                    FDoG_y=FDoG_y+self.score_list[n]*First_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
                    LoG_y=LoG_y+self.score_list[n]*LoG[sigma*times_of_sigma-x+n]
                    TDoG_y=TDoG_y+self.score_list[n]*Third_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
            elif x>len(self.score_list)-sigma*times_of_sigma-1:
                for n in range(x-sigma*times_of_sigma,len(self.score_list)):
                    y=y+self.score_list[n]*Gaussian[sigma*times_of_sigma-x+n]
                    FDoG_y=FDoG_y+self.score_list[n]*First_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
                    LoG_y=LoG_y+self.score_list[n]*LoG[sigma*times_of_sigma-x+n]
                    TDoG_y=TDoG_y+self.score_list[n]*Third_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
            return y,FDoG_y,LoG_y,TDoG_y

        def convolution_secondary(x,Gaussian,First_Derivative_of_Gaussian,LoG,Third_Derivative_of_Gaussian,sigma,times_of_sigma):    ### x is the index number of self.score_list, which is transmitted by parameter j in the following function "Gaussian_convolution_smoothing".
            #y=0
            #FDoG_y=0
            LoG_y=0
            #TDoG_y=0
            if x>=sigma*times_of_sigma and x<=len(self.score_list)-sigma*times_of_sigma-1:
                for n in range(x-sigma*times_of_sigma,x+sigma*times_of_sigma+1):
                    #y=y+self.score_list[n]*Gaussian[n-(x-sigma*times_of_sigma)]
                    #FDoG_y=FDoG_y+self.score_list[n]*First_Derivative_of_Gaussian[n-(x-sigma*times_of_sigma)]
                    LoG_y=LoG_y+self.score_list[n]*LoG[n-(x-sigma*times_of_sigma)]
                    #TDoG_y=TDoG_y+self.score_list[n]*Third_Derivative_of_Gaussian[n-(x-sigma*times_of_sigma)]
            elif x<sigma*times_of_sigma:
                for n in range(0,x+sigma*times_of_sigma+1):
                    #y=y+self.score_list[n]*Gaussian[sigma*times_of_sigma-x+n]
                    #FDoG_y=FDoG_y+self.score_list[n]*First_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
                    LoG_y=LoG_y+self.score_list[n]*LoG[sigma*times_of_sigma-x+n]
                    #TDoG_y=TDoG_y+self.score_list[n]*Third_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
            elif x>len(self.score_list)-sigma*times_of_sigma-1:
                for n in range(x-sigma*times_of_sigma,len(self.score_list)):
                    #y=y+self.score_list[n]*Gaussian[sigma*times_of_sigma-x+n]
                    #FDoG_y=FDoG_y+self.score_list[n]*First_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
                    LoG_y=LoG_y+self.score_list[n]*LoG[sigma*times_of_sigma-x+n]
                    #TDoG_y=TDoG_y+self.score_list[n]*Third_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
            return LoG_y  ##y,FDoG_y,LoG_y,TDoG_y

        ########
        Gaussian=self.ConvolutionParameters['Gaussian']
        First_Derivative_of_Gaussian=self.ConvolutionParameters['First_Derivative_of_Gaussian']
        LoG=self.ConvolutionParameters['LoG']
        Third_Derivative_of_Gaussian=self.ConvolutionParameters['Third_Derivative_of_Gaussian']
        sigma=self.ConvolutionParameters['sigma']
        ########
        times_of_sigma=self.ConvolutionParameters['times_of_sigma']
        ########
        secondary_Gaussian=self.ConvolutionParameters['secondary_Gaussian']
        secondary_First_Derivative_of_Gaussian=self.ConvolutionParameters['secondary_First_Derivative_of_Gaussian']
        secondary_LoG=self.ConvolutionParameters['secondary_LoG']
        secondary_Third_Derivative_of_Gaussian=self.ConvolutionParameters['secondary_Third_Derivative_of_Gaussian']
        secondary_sigma=self.ConvolutionParameters['secondary_sigma']
        ########
        for j in range(len(self.score_list)):
#            sys.stdout.write('\r        Coordinate:'+str(j*10+1) )
            Gaussian_result,FDoG_result,LoG_result,TDoG_result=convolution(j,Gaussian,First_Derivative_of_Gaussian,LoG,Third_Derivative_of_Gaussian,sigma,times_of_sigma)
            self.Gaussian_list.append(Gaussian_result)
            self.FDoG_list.append(FDoG_result)
            self.LoG_list.append(LoG_result)
            self.TDoG_list.append(TDoG_result)
            ########
            self.score_table.append([j*10+1,self.score_list[j],Gaussian_result,LoG_result,0])
            ########
            ##secondary_Gaussian_result,secondary_FDoG_result,secondary_LoG_result,secondary_TDoG_result=convolution(j,secondary_Gaussian,secondary_First_Derivative_of_Gaussian,secondary_LoG,secondary_Third_Derivative_of_Gaussian,secondary_sigma,times_of_sigma)
            secondary_LoG_result = convolution_secondary(j,secondary_Gaussian,secondary_First_Derivative_of_Gaussian,secondary_LoG,secondary_Third_Derivative_of_Gaussian,secondary_sigma,times_of_sigma)
            ##self.secondary_Gaussian_list.append(secondary_Gaussian_result)
            ##self.secondary_FDoG_list.append(secondary_FDoG_result)
            self.secondary_LoG_list.append(secondary_LoG_result)
            ##self.secondary_TDoG_list.append(secondary_TDoG_result)
            sys.stdout.flush()
        print('...... Convolution is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        print('\tChecking ==> Length of self.Gaussian_list:  ',len(self.Gaussian_list))##,'\tSecondary Length of self.Gaussian_list:  ',len(self.secondary_Gaussian_list))
        print('\tChecking ==> Length of self.FDoG_list:  ',len(self.FDoG_list))##,'\tSecondary Length of self.FDoG_list:  ',len(self.secondary_FDoG_list))
        print('\tChecking ==> Length of self.LoG_list:   ',len(self.LoG_list),'\tSecondary Length of self.LoG_list:   ',len(self.secondary_LoG_list))
        print('\tChecking ==> Length of self.TDoG_list:  ',len(self.TDoG_list))##,'\tSecondary Length of self.TDoG_list:  ',len(self.secondary_TDoG_list))




    def record_results(self):
        print('Record the detected nucleosomes in like-wiggle table ......')
        record_counting=0
        for j in range(len(self.final_nucleosome)):
            original_score_fragment=[]
            for nnp in range(self.final_nucleosome[j][2],self.final_nucleosome[j][3]+1):
                original_score_fragment.append(self.score_table[nnp][1])
            if len(original_score_fragment)==0:
                pass
            else:
                self.final_nucleosome[j][4]=max(original_score_fragment)
                self.final_nucleosome[j][5]=sum(original_score_fragment)
                for nnp in range(self.final_nucleosome[j][2],self.final_nucleosome[j][3]+1):
                    self.score_table[nnp][4]=max(original_score_fragment)    ### 在self.score_table里做记录
                record_counting=record_counting+1
        print('\tChecking ==> recorded nucleosomes:  ',record_counting)
        print('...... Nucleosomes are recorded in like-wiggle table.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))

         

    def write_wiggle_file(self):
        print('Generating like-wiggle file ......')
        wiggle=open(self.outputfile_wiggle,'w')
        wiggle.write('track type=like_wiggle_by_CHEN,Weizhong\nvariableStep chromosome='+self.chromosome+' span=10\n')
        ####wiggle.write('1:coordinate\t2:original_scoring\t3:Gaussian_convolution_smoothing\t4:First_Derivative_of_Gaussian_convolution\t5:LoG\t6:Tirst_Derivative_of_Gaussian_convolution\t7:Nucloesome_determination\n')
        wiggle.write('1:coordinate\t2:original_scoring\t3:Gaussian_convolution_smoothing\t4:LoG\t5:Adjusted_LoG\t6:Nucloesome_determination\n')
        line=0
        for k in range(len(self.score_list)):
            ####wiggle.write(str(k*10+1)+'\t'+str(self.score_list[k])+'\t'+str(self.Gaussian_list[k])+'\t'+str(self.FDoG_list[k])+'\t'+str(self.LoG_list[k])+'\t'+str(self.TDoG_list[k])+'\t'+str(self.secondary_LoG_list[k])+'\t'+str(score_file[k][4])+'\n')
            wiggle.write(str(k*10+1)+'\t'+str(self.score_list[k])+'\t'+str(self.Gaussian_list[k])+'\t'+str(self.LoG_list[k])+'\t'+str(self.secondary_LoG_list[k])+'\t'+str(self.score_table[k][4])+'\n')
            line=line+1
        wiggle.close()
        print('...... Wiggle file is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        print('\tChecking ==> Total line of the wiggle file:  ',line)














class NucleosomeAccuratePositioning(MainProgram):
    def __init__(self,FilesParameters,ConvolutionParameters,threshold):
        print(time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        self.inputfile = FilesParameters['inputfile_for_nucleosome_positioning']
        self.chromosome = FilesParameters['Chosen_Chromosome_Abbreviation']
        self.chrlength = FilesParameters['Chosen_Chromosome_Length']
        self.chrRecordingListLength = 0
        self.outputfile_wiggle = FilesParameters['outputfile_like_wig']
        self.outputfile_nucleosome = FilesParameters['outputfile_like_bed']
        print('Input file:  ',self.inputfile)
        ####print('Cell type:   ',self.celltype+' CD4+ T cells from human')
        print('Chromosome:  ',self.chromosome)
        print('Length of',self.chromosome,':  ',self.chrlength)
        print('Output wiggle file of nucleosome distribution: ',self.outputfile_wiggle)
        print('Output file of nucleosome collection:          ',self.outputfile_nucleosome)
        self.ConvolutionParameters = ConvolutionParameters
        self.threshold = threshold
        self.score_list = []
        self.tag_list = [ ]
        self.Gaussian_list = []
        self.FDoG_list = []
        self.LoG_list = []
        self.TDoG_list = []
        self.score_table = []    ### coordinate, self.score_list, self.Gaussian_list, self.LoG_list, column for recording
        ##self.secondary_Gaussian_list = []
        ##self.secondary_FDoG_list = []
        self.secondary_LoG_list = []
        ##self.secondary_TDoG_list = []
        print('======> ======> ======> ======> ======> ======>')
        #MainProgram.__init__(self,FilesParameters,ConvolutionParameters,threshold,chromosome_length)
        ######## 建立对象 self.inputfile
        ######## 建立对象 self.celltype
        ######## 建立对象 self.chromosome
        ######## 建立对象 self.chrlength
        ######## 建立对象 self.outputfile_path
        ######## 建立对象 self.outputfile_wiggle
        ######## 建立对象 self.outputfile_nucleosome
        ######## 建立对象 self.ConvolutionParameters=ConvolutionParameters
        ######## 建立对象 self.threshold=threshold
        ######## 建立对象 self.score_list=[]
        ######## 建立对象 self.Gaussian_list=[]
        ######## 建立对象 self.FDoG_list=[]
        ######## 建立对象 self.LoG_list=[]
        ######## 建立对象 self.TDoG_list=[]
        ######## 建立对象 self.score_table=[]
        ######## 建立对象 self.secondary_Gaussian_list=[]
        ######## 建立对象 self.secondary_FDoG_list=[]
        ######## 建立对象 self.secondary_LoG_list=[]
        ######## 建立对象 self.secondary_TDoG_list=[]

        if FilesParameters[ 'single_or_paired' ] == 's':
            MainProgram.score(self)
        elif FilesParameters[ 'single_or_paired' ] == 'p':
            MainProgram.score_for_Paired_end(self)
        
        MainProgram.Gaussian_convolution_smoothing(self)
        
      
        ####MainProgram.write_wiggle_file(self)
        ####MainProgram.write_nucleosome_file(self)

        ####def write_wiggle_file(self):
        print('Generating like-wiggle file ......')
        wiggle=open(self.outputfile_wiggle,'w')
        wiggle.write('track type=like_wiggle\nvariableStep chromosome='+self.chromosome+' span=10\n')
        ####wiggle.write('1:coordinate\t2:original_scoring\t3:Gaussian_convolution_smoothing\t4:First_Derivative_of_Gaussian_convolution\t5:LoG\t6:Tirst_Derivative_of_Gaussian_convolution\t7:Nucloesome_determination\n')
        wiggle.write('1:Coordinate\t2:Original_nucleosome_profile\t3:Gaussian_convolution_smoothing\t4:LoG\t5:Minor_LoG\t6:Tag_accumulation\t7:Nucloesomes\n')
        line=0
        for k in range(len(self.score_list)):
            ####wiggle.write(str(k*10+1)+'\t'+str(self.score_list[k])+'\t'+str(self.Gaussian_list[k])+'\t'+str(self.FDoG_list[k])+'\t'+str(self.LoG_list[k])+'\t'+str(self.TDoG_list[k])+'\t'+str(self.secondary_LoG_list[k])+'\t'+str(score_file[k][4])+'\n')
            wiggle.write(str(k*10+1)+'\t'+str(round(self.score_list[k],3))+'\t'+str(round(self.Gaussian_list[k],3))+'\t'+str(round(self.LoG_list[k],3))+'\t'+str(round(self.secondary_LoG_list[k],3))+'\t'+str(self.tag_list[k])+'\t'+str(round(self.score_table[k][4],3))+'\n')
            line=line+1
        wiggle.close()
        print('...... Wiggle file is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        print('\tChecking ==> Total line of the wiggle file:  ',line)

       
        print('...... Nucleosome collection file is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        print('\tChecking ==> Total line of the collection file:  ',line)


        print('\n')



        





def count_chromosome_length( input_file ):
    print('Determine chromosome length ......')
    max_coordinate = 0
    inputTEXT = open( input_file , 'r' )
    while True:
        try:
            line2list = next( inputTEXT ).split('\n')[0].split('\r')[0].split('\t')
            if line2list[ 5 ] == '+':
                right_coordinate = int( line2list[1] ) + 150
            elif line2list[ 5 ] == '-':
                right_coordinate = int( line2list[2] )
            ########
            if right_coordinate > max_coordinate:
                max_coordinate = right_coordinate
            else:
                pass
        except StopIteration:
            break
    inputTEXT.close()
    print('...... Default chromosome length is determined.\n')
    ##return max_coordinate+100
    return max_coordinate









class input_data_pretreatment_for_single_chromosome:
    def __init__( self , parameters_pretreatment_for_single_chromosome ):
        self.chr_recording_dict = { }
        self.chr_recording_list = [ ]
        self.chromosomeName     = parameters_pretreatment_for_single_chromosome[ 'chromosome' ]
        self.chromosome_lines         = 0
        self.chromosome_maxCoordinate = 0
        self.chromosome_lengthSetting = 0
        self.OutputFile_Tags    = parameters_pretreatment_for_single_chromosome[ 'Output_path' ] + self.chromosomeName + '.bed'
        self.OutputFile_summary = parameters_pretreatment_for_single_chromosome[ 'OutputFile_summary' ]
        ########
        if parameters_pretreatment_for_single_chromosome['single_or_paired_setting'] == 's':
            self.do_pretreatment( parameters_pretreatment_for_single_chromosome )
        elif parameters_pretreatment_for_single_chromosome['single_or_paired_setting'] == 'p':
            self.do_pretreatment_for_paired_end( parameters_pretreatment_for_single_chromosome )
        ########
        self.write_down_summary( parameters_pretreatment_for_single_chromosome )

    def do_pretreatment( self , parameters_pretreatment_for_single_chromosome ):
        count_raw_line = 0 #### 计数
        print('Read raw data: ' , parameters_pretreatment_for_single_chromosome[ 'Raw_input_file' ] )
        print('File for preliminary data: ' , self.OutputFile_Tags )
        inputTEXT  = open( parameters_pretreatment_for_single_chromosome[ 'Raw_input_file' ] , 'r' ) #### 只读
        outputTEXT = open( self.OutputFile_Tags , 'w' )
        while True:
            try:
                original_line = next( inputTEXT )
                count_raw_line = count_raw_line + 1  #### 计数
 #               sys.stdout.write('\r    Raw line: '+str(count_raw_line) )
                line2list = original_line.split('\n')[0].split('\r')[0].split('\t')
                ####
                if line2list[ 5 ] == '+':
                    right_coordinate = int( line2list[1] ) + 150
                elif line2list[ 5 ] == '-':
                    right_coordinate = int( line2list[2] )
                ####
                chromosome = line2list[ 0 ]
                self.chr_recording_dict[ chromosome ] = ''
                ##
                if self.chromosomeName == chromosome:
                    outputTEXT.write( original_line )
                    self.chromosome_lines = self.chromosome_lines + 1
                    if right_coordinate > self.chromosome_maxCoordinate:
                        self.chromosome_maxCoordinate = right_coordinate
                    else:
                        pass
                else:
                    pass
                ####
                sys.stdout.flush()
            except StopIteration:
                break
        inputTEXT.close() #### 读毕
        outputTEXT.close()
        print('\tChecking ==> Lines of raw data: ' , count_raw_line )
        print('\tChecking ==> Lines for ' , self.chromosomeName , ': ' , self.chromosome_lines )
        if parameters_pretreatment_for_single_chromosome[ 'length' ] == 'none':
            self.chromosome_lengthSetting = self.chromosome_maxCoordinate
        elif parameters_pretreatment_for_single_chromosome[ 'length' ] != 'none':
            self.chromosome_lengthSetting = parameters_pretreatment_for_single_chromosome[ 'length' ]
        ####
        self.chr_recording_list = list( self.chr_recording_dict.keys() )

    def do_pretreatment_for_paired_end( self , parameters_pretreatment_for_single_chromosome ):
        count_raw_line = 0 #### 计数
        print('Read raw data: ' , parameters_pretreatment_for_single_chromosome[ 'Raw_input_file' ] )
        print('File for preliminary data: ' , self.OutputFile_Tags )
        inputTEXT  = open( parameters_pretreatment_for_single_chromosome[ 'Raw_input_file' ] , 'r' ) #### 只读
        outputTEXT = open( self.OutputFile_Tags , 'w' )
        while True:
            try:
                original_line = next( inputTEXT )
                count_raw_line = count_raw_line + 1  #### 计数
#                sys.stdout.write('\r    Raw line: '+str(count_raw_line) )
                line2list = original_line.split('\n')[0].split('\r')[0].split('\t')
                ####
                ##if line2list[ 5 ] == '+':
                ##    right_coordinate = int( line2list[1] ) + 150
                ##elif line2list[ 5 ] == '-':
                ##    right_coordinate = int( line2list[2] )
                right_coordinate = int( line2list[2] )
                ####
                chromosome = line2list[ 0 ]
                self.chr_recording_dict[ chromosome ] = ''
                ##
                if self.chromosomeName == chromosome:
                    outputTEXT.write( original_line )
                    self.chromosome_lines = self.chromosome_lines + 1
                    if right_coordinate > self.chromosome_maxCoordinate:
                        self.chromosome_maxCoordinate = right_coordinate
                    else:
                        pass
                else:
                    pass
                ####
                sys.stdout.flush()
            except StopIteration:
                break
        inputTEXT.close() #### 读毕
        outputTEXT.close()
        print('\tChecking ==> Lines of raw data: ' , count_raw_line )
        print('\tChecking ==> Lines for ' , self.chromosomeName , ': ' , self.chromosome_lines )
        if parameters_pretreatment_for_single_chromosome[ 'length' ] == 'none':
            self.chromosome_lengthSetting = self.chromosome_maxCoordinate
        elif parameters_pretreatment_for_single_chromosome[ 'length' ] != 'none':
            self.chromosome_lengthSetting = parameters_pretreatment_for_single_chromosome[ 'length' ]
        ####
        self.chr_recording_list = list( self.chr_recording_dict.keys() )

    def write_down_summary( self , parameters_pretreatment_for_single_chromosome ):
        print('Write down summary: ' , self.OutputFile_summary )
        outputTEXT = open( self.OutputFile_summary , 'w' )
        outputTEXT.write( ('\t').join( [ 'Chromosome' , 'Tags' , 'Max_coordinate' , 'Recorded_chromosome_length' ] ) + '\n' )
        outputTEXT.write( ('\t').join( [ self.chromosomeName , str(self.chromosome_lines ) , str(self.chromosome_maxCoordinate) , str(self.chromosome_lengthSetting) ] ) + '\n' )
        outputTEXT.close()








class input_data_pretreatment:
    def __init__( self , parameters_pretreatment ):
        self.chr_list = [ ]
        self.chr2maxCoordinate_dict = { }
        self.chr2lengthSetting_dict = { }
        self.chr2lines_dict = { }
        self.chr2outputfile_dict = { }
        ########
        if parameters_pretreatment['single_or_paired_setting'] == 's':
            self.do_pretreatment( parameters_pretreatment )
        elif parameters_pretreatment['single_or_paired_setting'] == 'p':
            self.do_pretreatment_for_paired_end( parameters_pretreatment )
        ########
        self.write_down_summary( parameters_pretreatment )

    def do_pretreatment( self , parameters_pretreatment ):
        chr2outputTEXT_dict = { }
        count_raw_line = 0 #### 计数
        print('Read raw data: ' , parameters_pretreatment['Raw_input_file'] )
        print('Path to record preliminary data: ' , parameters_pretreatment['Output_path'] )
        inputTEXT = open( parameters_pretreatment['Raw_input_file'] , 'r' ) #### 只读
        while True:
            try:
                original_line = next( inputTEXT )
                count_raw_line = count_raw_line + 1  #### 计数
#                sys.stdout.write('\r    Raw line: '+str(count_raw_line) )
                line2list = original_line.split('\n')[0].split('\r')[0].split('\t')
                ####
                if line2list[ 5 ] == '+':
                    right_coordinate = int( line2list[1] ) + 150
                elif line2list[ 5 ] == '-':
                    right_coordinate = int( line2list[2] )
                ####
                chromosome = line2list[ 0 ]
                if chromosome not in self.chr2outputfile_dict.keys():
                    self.chr2outputfile_dict[ chromosome ] = parameters_pretreatment['Output_path']+chromosome+'.bed'
                    chr2outputTEXT_dict[ chromosome ] = open( self.chr2outputfile_dict[ chromosome ] , 'w' )
                    self.chr2maxCoordinate_dict[ chromosome ] = 0
                    self.chr2lines_dict[ chromosome ] = 0
                    self.chr_list.append( chromosome )
                elif chromosome in self.chr2outputfile_dict.keys():
                    pass
                ################
                chr2outputTEXT_dict[ chromosome ].write( original_line )
                self.chr2lines_dict[ chromosome ] = self.chr2lines_dict[ chromosome ] + 1
                ####
                if right_coordinate > self.chr2maxCoordinate_dict[ chromosome ]:
                    self.chr2maxCoordinate_dict[ chromosome ] = right_coordinate
                else:
                    pass
                ####
                sys.stdout.flush()
            except StopIteration:
                break
        inputTEXT.close() #### 读毕
        print('\tChecking ==> Lines of raw data: ' , count_raw_line )
        ########
        for chromosome in self.chr_list:
            chr2outputTEXT_dict[ chromosome ].close()
            ##self.chr2lengthSetting_dict[ chromosome ] = self.chr2maxCoordinate_dict[ chromosome ] + 100
            self.chr2lengthSetting_dict[ chromosome ] = self.chr2maxCoordinate_dict[ chromosome ]
        ########

    def do_pretreatment_for_paired_end( self , parameters_pretreatment ):
        chr2outputTEXT_dict = { }
        count_raw_line = 0 #### 计数
        print('Read raw data: ' , parameters_pretreatment['Raw_input_file'] )
        print('Path to record preliminary data: ' , parameters_pretreatment['Output_path'] )
        inputTEXT = open( parameters_pretreatment['Raw_input_file'] , 'r' ) #### 只读
        while True:
            try:
                original_line = next( inputTEXT )
                count_raw_line = count_raw_line + 1  #### 计数
#                sys.stdout.write('\r    Raw line: '+str(count_raw_line) )
                line2list = original_line.split('\n')[0].split('\r')[0].split('\t')
                ####
                ##if line2list[ 5 ] == '+':
                ##    right_coordinate = int( line2list[1] ) + 150
                ##elif line2list[ 5 ] == '-':
                ##    right_coordinate = int( line2list[2] )
                right_coordinate = int( line2list[2] )
                ####
                chromosome = line2list[ 0 ]
                if chromosome not in self.chr2outputfile_dict.keys():
                    self.chr2outputfile_dict[ chromosome ] = parameters_pretreatment['Output_path']+chromosome+'.bed'
                    chr2outputTEXT_dict[ chromosome ] = open( self.chr2outputfile_dict[ chromosome ] , 'w' )
                    self.chr2maxCoordinate_dict[ chromosome ] = 0
                    self.chr2lines_dict[ chromosome ] = 0
                    self.chr_list.append( chromosome )
                elif chromosome in self.chr2outputfile_dict.keys():
                    pass
                ################
                chr2outputTEXT_dict[ chromosome ].write( original_line )
                self.chr2lines_dict[ chromosome ] = self.chr2lines_dict[ chromosome ] + 1
                ####
                if right_coordinate > self.chr2maxCoordinate_dict[ chromosome ]:
                    self.chr2maxCoordinate_dict[ chromosome ] = right_coordinate
                else:
                    pass
                ####
                sys.stdout.flush()
            except StopIteration:
                break
        inputTEXT.close() #### 读毕
        print('\tChecking ==> Lines of raw data: ' , count_raw_line )
        ########
        for chromosome in self.chr_list:
            chr2outputTEXT_dict[ chromosome ].close()
            ##self.chr2lengthSetting_dict[ chromosome ] = self.chr2maxCoordinate_dict[ chromosome ] + 100
            self.chr2lengthSetting_dict[ chromosome ] = self.chr2maxCoordinate_dict[ chromosome ]
        ########

    def write_down_summary( self , parameters_pretreatment ):
        print('Write down summary: ' , parameters_pretreatment['OutputFile_summary'] )
        outputTEXT = open( parameters_pretreatment['OutputFile_summary'] , 'w' )
        outputTEXT.write( ('\t').join( [ 'Chromosome' , 'Tags' , 'Max_coordinate' , 'Recorded_chromosome_length' ] ) + '\n' )
        for chromosome in sorted( self.chr_list ):
            outputTEXT.write( ('\t').join( [ chromosome , str(self.chr2lines_dict[ chromosome ]) , str(self.chr2maxCoordinate_dict[ chromosome ]) , str(self.chr2lengthSetting_dict[ chromosome ]) ] ) + '\n' )
        outputTEXT.close()











def Merge_overall_results( parameters_Merge_results ):
    print('Gather all nucleosome detection results: ' , parameters_Merge_results['OutputFile'] )
    outputTEXT = open( parameters_Merge_results['OutputFile'] , 'w' )
    total_nucleosome____dict = { }
    headline2 = ''
    for chromosome in sorted( parameters_Merge_results['Chromosome_list'] ):
        InputFile = parameters_Merge_results['Input_Files'] + '_' + chromosome + '.like_bed'
        inputTEXT = open( InputFile , 'r' ) #### 只读
        headline1 = next( inputTEXT ) #### headline1
        outputTEXT.write( headline1 ) #### 记录
        headline2 = next( inputTEXT ) #### headline2
        inputTEXT.close() #### 读毕
    outputTEXT.write( headline2 )     #### 记录
    ########
    for chromosome in sorted( parameters_Merge_results['Chromosome_list'] ):
        InputFile = parameters_Merge_results['Input_Files'] + '_' + chromosome + '.like_bed'
        print('Read from: ' , InputFile )
        inputTEXT = open( InputFile , 'r' ) #### 只读
        next( inputTEXT ) #### headline1
        next( inputTEXT ) #### headline2
        while True:
            try:
                original_line = next( inputTEXT )
                outputTEXT.write( original_line )     #### 记录
            except StopIteration:
                break
        inputTEXT.close() #### 读毕
    outputTEXT.close()
    print('...... All the nucleosome results are gathered.')
    print('\n')











if __name__=='__main__':
    
    print('\n\nProgram to detect nucleosomes from MNase-seq data.')
    print('\nThe program should run in Python3.\n')
    usage='usage: iNPS.py [options]'
    parser = OptionParser( usage = '%prog  -i /path/inputfile  -o /path/outputfile  -c chromosome_name  -l chromosome_length  --s_p single_or_paired_end_data' , version = '%prog Version:1.2.2' )
    parser.add_option('-i', '--input',      dest='input_file',                        type='string', help='"/path/filename"  INPUT_FILE file of single-end sequencing tags in a standard BED format ( chromosome <tab> start <tab> end <tab> name <tab> score <tab> strand ), or paired-end tags in a 3-column BED format ( chromosome <tab> start <tab> end ).')
    parser.add_option('-o', '--output',     dest='output_file',                       type='string', help='"/path/filename"  Here, the name extension is unnecessary. Software will output two result files, "filename_[ChromosomeName].like_bed" and "filename_[ChromosomeName].like_wig", to record coordinates and profiles of detected nucleosomes respectively. The chromosome name will be added as suffix in the file names. If your detect nucleosomes on multiple chromosomes, for each chromosome, software will output two result files "filename_[ChromosomeName].like_bed" and "filename_[ChromosomeName].like_wig" respectively. And finally, a file "filename_Gathering.like_bed" will gather the detected nucleosomes on every chromosome. Note that a path "/path/filename/" or "/path/filename_[ChromosomeName]/" will be built to record the preliminary and intermediate data.')
    parser.add_option('-c', '--chrname',    dest='chromosome_name',                   type='string', help='Specify the name (or abbreviation) of the chromosome, if you would like to do nucleosome detection ONLY on ONE single chromosome. For nucleosome detection on multiple chromosomes, please do NOT use this parameter. That is, if your do NOT use this parameter, software will detect nucleosome on each chromosome ONE-BY-ONE in the input data as default.')
    parser.add_option('-l', '--chrlength',  dest='chromosome_length',                 type='int',    help='The length of the chromosome specified by parameter "-c" or "--chrname". ONLY used for nucleosome detection on ONE single chromosome (parameter "-c" or "--chrname" is setted). If you do NOT use this parameter, software will find the maximum coordinate in the input data to represent the chromosome length as default. For nucleosome detection on multiple chromosomes, please do NOT use this parameter. The length of each chromosome will be determined by the tag with maximum coordinate of the corresponding chromosome respectively.')
    parser.add_option('--s_p',              dest='single_or_paired_end',              type='string', default='s' , help='"s" or "p". [Default = s] Set to "p" if the input data is paired-end tags. Otherwise, set to "s" or use the default setting if the input data is single-end tags.')
    parser.add_option('--pe_max',           dest='superior_limit_of_paired_end_tags', type='int',    default=200 , help='The superior limit of the length of paired-end tags. [Default = 200] The tags longer than the cutoff will be ignored. This parameter is ONLY available for paired-end sequencing data. Please avoid using too large value.')
    parser.add_option('--pe_min',           dest='inferior_limit_of_paired_end_tags', type='int',    default=100 , help='The inferior limit of the length of paired-end tags. [Default = 100] The tags shorter than the cutoff will be ignored. This parameter is ONLY available for paired-end sequencing data. Please avoid using too small value.')
    (options, args) = parser.parse_args( )
    ########
    if options.input_file == None:
        parser.error('-h for help or provide the input file name!')
    else:
        pass
    ########
    if options.output_file == None:
        parser.error('-h for help or provide the output file name!')
    else:
        pass
    ########
    if options.single_or_paired_end == 's':
        FilesParameters[ 'single_or_paired' ] = 's'
    elif options.single_or_paired_end == 'p':
        FilesParameters[ 'single_or_paired' ] = 'p'
        FilesParameters[ 'pe_max' ] = options.superior_limit_of_paired_end_tags
        FilesParameters[ 'pe_min' ] = options.inferior_limit_of_paired_end_tags
    else:
        parser.error('-h for help or correctly specify the type of input data: single-end ("s") or paired_end ("p")')





    ########
    print('Prepare Gaussian function and Laplacian of Gaussian operator for convolution ......')
    ConvolutionParameters=Gaussian_profile(PreliminaryConvolutionParameters)
    ####print('\tDiscrete Gaussian function:\t(Range:',len(ConvolutionParameters['Gaussian']),')')
    ####print('\t\t','Gaussian','\t\t','1st_Der','\t\t','2nd_Der','\t\t','3rd_Der')
    ####for i in range(len(ConvolutionParameters['Gaussian'])):
    ####    print('\t\t',round(ConvolutionParameters['Gaussian'][i],6),'\t\t',round(ConvolutionParameters['First_Derivative_of_Gaussian'][i],6),'\t\t',round(ConvolutionParameters['LoG'][i],6),'\t\t',round(ConvolutionParameters['Third_Derivative_of_Gaussian'][i],6))
    ####print('\tDiscrete secondary Gaussian function:\t(Range:',len(ConvolutionParameters['secondary_Gaussian']),')')
    ####print('\t\t','Gaussian','\t\t','1st_Der','\t\t','2nd_Der','\t\t','3rd_Der')
    ####for i in range(len(ConvolutionParameters['secondary_Gaussian'])):
    ####    print('\t\t',round(ConvolutionParameters['secondary_Gaussian'][i],6),'\t\t',round(ConvolutionParameters['secondary_First_Derivative_of_Gaussian'][i],6),'\t\t',round(ConvolutionParameters['secondary_LoG'][i],6),'\t\t',round(ConvolutionParameters['secondary_Third_Derivative_of_Gaussian'][i],6))
    print('...... The preparation of Gaussian convolution parameters is finished.') 
    print('\n')
    ########





    if options.chromosome_name != None:
        print('Do nucleosome detection on chromosome: ' , options.chromosome_name )
        print('\n')
        ######## inter-mediate path ########
        Intermediate_Path = options.output_file+'_'+options.chromosome_name+'/' #### inter-mediate
        if os.path.exists( Intermediate_Path ):
            pass
        else:
            os.makedirs( Intermediate_Path )
        ######## Data pretreatment ########
        print('Do data pretreatment ......')
        parameters_pretreatment_for_single_chromosome = { 'Raw_input_file': options.input_file,
                                                          'Output_path':    Intermediate_Path,
                                                          'OutputFile_summary': Intermediate_Path+'InputData_Summary.txt',
                                                          'chromosome':   options.chromosome_name ,
                                                          'single_or_paired_setting': options.single_or_paired_end,
                                                          }
        if options.chromosome_length == None:
            parameters_pretreatment_for_single_chromosome[ 'length' ] = 'none'
        elif options.chromosome_length != None:
            parameters_pretreatment_for_single_chromosome[ 'length' ] = options.chromosome_length
        InputData_Summary_for_single_chromosome = input_data_pretreatment_for_single_chromosome( parameters_pretreatment_for_single_chromosome ) #### get the summary of input data
        ####    InputData_Summary_for_single_chromosome.chr_recording_dict
        ####    InputData_Summary_for_single_chromosome.chr_recording_list
        ####    InputData_Summary_for_single_chromosome.chromosomeName
        ####    InputData_Summary_for_single_chromosome.chromosome_lines
        ####    InputData_Summary_for_single_chromosome.chromosome_maxCoordinate
        ####    InputData_Summary_for_single_chromosome.chromosome_lengthSetting
        ####    InputData_Summary_for_single_chromosome.OutputFile_Tags
        ####    InputData_Summary_for_single_chromosome.OutputFile_summary
        print('...... Data pretreatment finished.')
        print('\n')
        ########
        if ( options.chromosome_name in InputData_Summary_for_single_chromosome.chr_recording_list ) and ( InputData_Summary_for_single_chromosome.chromosome_lines > 0 ):
            FilesParameters[ 'Chosen_Chromosome_Abbreviation' ] = options.chromosome_name
            FilesParameters[ 'inputfile_for_nucleosome_positioning' ] = InputData_Summary_for_single_chromosome.OutputFile_Tags
            FilesParameters[ 'outputfile_like_bed' ] = options.output_file + '_' + options.chromosome_name + '.like_bed'
            FilesParameters[ 'outputfile_like_wig' ] = options.output_file + '_' + options.chromosome_name + '.like_wig'
            if options.chromosome_length == None:
                FilesParameters[ 'Chosen_Chromosome_Length' ] = InputData_Summary_for_single_chromosome.chromosome_lengthSetting
            elif options.chromosome_length != None:
                FilesParameters[ 'Chosen_Chromosome_Length' ] = options.chromosome_length
            #### 开始执行函数 ####
            #### 开始执行函数 ####
            #### 开始执行函数 ####
            NucleosomeAccuratePositioning( FilesParameters , ConvolutionParameters , threshold )
            #### 执行函数完毕 ####
            #### 执行函数完毕 ####
            #### 执行函数完毕 ####
        else:
            print('Please provide correct chromosome name, or choose chromosome among: ' , (',').join(sorted(InputData_Summary_for_single_chromosome.chr_recording_list)) )
    ################
    ################
    elif options.chromosome_name == None:
        ######## inter-mediate path ########
        Intermediate_Path = options.output_file+'/' #### inter-mediate
        if os.path.exists( Intermediate_Path ):
            pass
        else:
            os.makedirs( Intermediate_Path )
        ######## Data pretreatment ########
        print('Do data pretreatment ......')
        parameters_pretreatment = { 'Raw_input_file': options.input_file,
                                    'Output_path':    Intermediate_Path,
                                    'OutputFile_summary': Intermediate_Path+'InputData_Summary.txt',
                                    'single_or_paired_setting': options.single_or_paired_end,
                                    }
        InputData_Summary = input_data_pretreatment( parameters_pretreatment ) #### get the summary of input data
        ####    InputData_Summary.chr_list = [ ]
        ####    InputData_Summary.chr2maxCoordinate_dict = { }
        ####    InputData_Summary.chr2lengthSetting_dict = { }
        ####    InputData_Summary.chr2lines_dict = { }
        ####    InputData_Summary.chr2outputfile_dict = { }
        print('...... Data pretreatment finished.')
        print('\n')
        ########
        for each_chromosome in sorted(InputData_Summary.chr_list):
            FilesParameters[ 'Chosen_Chromosome_Abbreviation' ] = each_chromosome
            FilesParameters[ 'inputfile_for_nucleosome_positioning' ] = InputData_Summary.chr2outputfile_dict[ each_chromosome ]
            FilesParameters[ 'outputfile_like_bed' ] = options.output_file + '_' + each_chromosome + '.like_bed'
            FilesParameters[ 'outputfile_like_wig' ] = options.output_file + '_' + each_chromosome + '.like_wig'
            FilesParameters[ 'Chosen_Chromosome_Length' ] = InputData_Summary.chr2lengthSetting_dict[ each_chromosome ]
            #### 开始执行函数 ####
            #### 开始执行函数 ####
            #### 开始执行函数 ####
            NucleosomeAccuratePositioning( FilesParameters , ConvolutionParameters , threshold )
            #### 执行函数完毕 ####
            #### 执行函数完毕 ####
            #### 执行函数完毕 ####
        ######## 合并 nucleosomes ########
        parameters_Merge_results = { 'Input_Files': options.output_file,
                                     'Chromosome_list': InputData_Summary.chr_list,
                                     'OutputFile':  options.output_file + '_Gathering.like_bed',
                                     }
        print('\n')
        Merge_overall_results( parameters_Merge_results )
    ########





print('\nThe end.\n\n')



