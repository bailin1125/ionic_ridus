##这个程序的目的是输入一个离子半径表，输出一个填充之后的离子半径表
import numpy as np
from scipy.optimize import leastsq
import os



#定义元素符号和文件名
atom_name=["D","H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg"]
file_name="离子半径表_clean.prn"
out_name="离子半径表_full_fill"
iteration="Test the number of iteration"
##用来测试的数据
Xi=np.array([8.19,2.72,6.39,8.71,4.7,2.66,3.78])
Yi=np.array([7.01,2.78,6.47,6.71,4.1,4.23,4.05])
class ionic_ridus:
    "离子半径的类定义"
    def __init__(self,str_name):       
        self.atom_num=atom_name.index(str_name[0])    ## ionic atomic number
        self.val_state=int(str_name[1]) ## ionic valence state    
        self.coor_num=int(str_name[2])    ## ionic coordination number 
        self.structure_type=int(str_name[3])    ## ionic structure type, 0=common, 1=square, 2=pyrometric cone
        self.spin_stat=int(str_name[4])   ## ionic spin state, 0=common, 1=high spin, 2=low spin
        self.radius=float(str_name[5])     ## ionic radius
    def print(self):
        print("the atom_num={}".format(self.atom_num))
        print("the coor_num={}".format(self.coor_num))
        print("the val_state={}".format(self.val_state))
        print("the radius={}".format(self.radius))
        print("the structurs type ={}".format(self.structure_type))
        print("the spin_state={}".format(self.spin_stat))


class ionic_plus:
    ##用来扩展的离子半径   
    type_num=0
    def __init__(self,atom_num,valence):
        self.atom_num=atom_num
        self.valence=valence
        self.coor=[]
        self.ridus=[]
        self.org_num=0##填充前储存的键长信息
        self.change=False
        ionic_plus.type_num+=1
    def puls_info(self,coor,ridus):
        self.coor.append(float(coor))
        self.ridus.append(float(ridus))
        self.org_num+=1
    def fll_fill(self):##通过最小二乘法补充完配位数对应的半径信息
        if(len(self.coor)<2):
            return
        result=get_least_logirm(error,p0,arggs=(np.array(self.coor),np.array(self.ridus)))
        k,b=result[0]
        min_temp=int(min(self.coor))
        max_temp=int(max(self.coor))
        for i in range(min_temp,max_temp):
            try: 
                self.coor.index(i)
            except:
                self.change=True
                self.coor.append(i)
                self.ridus.append(k*i+b)   
    def print_info(self):
        print("for {} element, val is {}\n".format(atom_name[self.atom_num],self.valence))
        print("{}\n".format(self.coor))
        print("{}\n".format(self.ridus))
  
##最小二乘法的相关函数
def func(p,x):
    k,b=p
    return k*x+b

def error(p,x,y):
    #print(s)
    return func(p,x)-y

def get_least_logirm(error,p0,arggs):
    para=leastsq(error,p0,arggs)
    k,b=para[0]
    #print("k={},b={}".format(k,b))
    return para

def output_result(ridus,ionic_plus_real,file_name):
    out=open(file_name,"w")
    for i in range(len(ridus)):
        out.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{:.2f}\n".format(atom_name[ridus[i].atom_num],ridus[i].val_state,
                  ridus[i].coor_num,ridus[i].structure_type,ridus[i].spin_stat,
                  ridus[i].radius))

    out.close()
    
    ##然后输出的是修正后的结果
    out=open(file_name+"_corect","w")
    for plus in ionic_plus_real:
        if(plus.change is True):
            for  i in range(plus.org_num,len(plus.coor)):                
                out.write("{}\t\t{}\t\t{}\t\t{:.2f}\n".format(atom_name[plus.atom_num],plus.valence,
                 int(plus.coor[i]),plus.ridus[i]))
    out.close()

    
##读取文件
def read_ridus(filename):
    ridus=[]
    i=0
    fo=open(filename,"r")
    for line in fo:
        read=line.split()        
        ridus.append(ionic_ridus(read)) 
        #print("\n")
        #ridus[i].print()
        i+=1
    fo.close()
    return ridus
def sort_key(ionic):
    return ionic.atom_num
if __name__=="__main__":
    ##首先读取文件中的相关信息
    ridus=read_ridus(file_name)
    ##先看一下那样写有没有问题
    for i in range(len(ridus)-1):
        if(ridus[i].atom_num==ridus[i+1].atom_num and ridus[i].val_state==ridus[i+1].val_state and 
           ridus[i].coor_num==ridus[i+1].coor_num and ( ridus[i].structure_type!=ridus[i+1].structure_type or ridus[i].spin_stat!=ridus[i+1].spin_stat)):
            print("the element {} has coincendenc!".format(atom_name[ridus[i].atom_num]))
        
    print("test completed!")
    ##简单的测试一下最小二乘法的威力
    p0=[1,2]
    arggs=(Xi,Yi)
    result=get_least_logirm(error,p0,arggs)


    ##然后开始进行最小二乘法扩充
    print("开始填充ionic_plus对象")
    ionic_plus_real=[]
    flag=[]##标记每个原始离子半径是否被用到了
    for i in range(1,len(ridus)+1):
        if(i==len(ridus)):
            break
        now_atom=(ridus[i].atom_num,ridus[i].val_state,ridus[i].coor_num,ridus[i].radius)
        coor=ridus[i].coor_num
        riduss=ridus[i].radius
        if(ridus[i].spin_stat !=0):
            continue
        if(ridus[i].atom_num!= ridus[i-1].atom_num  or ridus[i].val_state!= ridus[i-1].val_state):
            ionic_plus_real.append(ionic_plus(now_atom[0],now_atom[1]))
            ionic_plus_real[-1].puls_info(now_atom[2],now_atom[3])                     
        else:
            ionic_plus_real[ionic_plus.type_num-1].puls_info(coor,riduss)
    ##开始进行最小二乘法
    for i in range(len(ionic_plus_real)):##开始真正填充每个对象了
        #ionic_plus_real[i].print_info()
        if(len(ionic_plus_real[i].coor)>1):
            ionic_plus_real[i].fll_fill()            
            #ionic_plus_real[i].print_info()
        #os.system("pause")
    ridus.sort(key=sort_key)
    ionic_plus_real.sort(key=sort_key)

    #填充完对象开始写文件了
    print("开始写入文件！")
    output_result(ridus,ionic_plus_real,out_name)
    print("写入文件完毕，请检查路径中的文件")

