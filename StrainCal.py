# pylint: disable=C0111,R0902,R0902,W0611,W0301,C0103,C0411
from ase.io import read;
from ase.visualize import view;
from ase import Atoms
from scipy.interpolate import CubicSpline,interp1d,PchipInterpolator,Akima1DInterpolator,KroghInterpolator
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
import numpy as num
import numpy as np;
from ase.io import read,write;
import os;
from configparser import ConfigParser;
import sys
import math
import pandas as pd


class QERun:
    def __init__(self,inputName = "test.in",runCommand = ["wsl","mpirun","-np","4","pw.x"]):
        self.file = open(inputName);
        self.nat = 0;
        self.coords = list();
        self.atomsLabels = list();
        self.cell = list();
        self.unit = "crystal";
        self.cellUnit = "angstrom";
        self.text = "";
        self.atoms = Atoms();
        self.atoms = read(inputName);
        self.currentValue = 0;
        self.inputText = "";
        self.runCommand = runCommand;
        self.currentCell = [];
        self.type = 'all'
        self.config = ConfigParser();
        self.config.read("config.ini");
        self.mainVol = 0;
        self.properties = pd.DataFrame(columns = ["volume","deltaV","alpha","beta","gamma","energy"])
        self.xy = pd.DataFrame(columns = ['x','y'])
        

    def volumeCalc(self,cell):
        x = cell;
        s = np.cross(x[0],x[1]);
        v = np.dot(s,x[2]);
        return v;

    def readFile(self):
        readCoords = False;
        lineCounter = 1;
        readCell = False;
        text = ""
        for line in self.file.readlines():
            words = line.replace(" ","").replace(",","").split("=");
            if "nat" in words[0]:
                self.nat = int(words[1]);
            if line.find("ATOMIC_POSITIONS") >-1:
                self.unit = line.split()[-1];
                readCoords = True;
                lineCounter = 1;
                continue;
            if readCoords and lineCounter <= self.nat:
                w = line.split();
                self.coords.append([float(w[1]), float(w[2]), float(w[3])]);
                self.atomsLabels.append(w[0])

                lineCounter += 1;
                continue;
            if line.find("CELL_PARAMETERS") > -1:
                w = line.replace(",","").split();
                self.cellUnit = w[-1];
                readCell = True;
                lineCounter = 1;
                continue;
            if readCell and lineCounter <=3:
                self.cell.append([float(w) for w in line.split()]);
                lineCounter +=1;
                continue;
            text += line;
        self.text = text;
        self.file.close()
        self.mainVol = self.volumeCalc(self.cell)
        #view(self.atoms)

    def ase2Cor(self):
            if (self.unit == "crystal"):
                self.coords = self.atoms.get_scaled_positions();
            if (self.unit == "angstrom"):
                self.coords = self.atoms.get_positions();
    def changeCoorCells(self,newCell):
        fracPos = self.atoms.get_scaled_positions();
        atoms = self.atoms;
        atoms.set_cell(newCell);
        atoms.set_scaled_positions(fracPos);
        return atoms.get_scaled_positions();

    def rotation_matrix(self,vector, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.dot(vector, np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]]))

    def cellMaker(self,maximum = 0.1,steps = 0.01,direction="all"):
        print(np.arange(-maximum,maximum,steps));
        for step in np.arange(-maximum,maximum+steps,steps):
            step = np.round(step,self.config.getint("run","rounding_digit_step"));
            self.currentValue = step;
            self.type = direction;
            tempCell = [[[] for i in range(3)] for i in range(3)];#Creating an empty list
            a = self.cell[0]
            b = self.cell[1]
            c = self.cell[2]
            for i in range(0,3):
                for j in range(0,3):
                    tempCell[i][j] = self.cell[i][j];
            if (direction=='all'):
                for i in range(0,3):
                    for j in range (0,3):
                        tempCell[i][j] = self.cell[i][j]+self.cell[i][j]*step;
            if (direction=='x'):
                for i in range(0,3):
                    tempCell[i][0] = self.cell[i][0] + self.cell[i][0]*step;
            if (direction=='xy'):
                for i in range(0,3):
                    tempCell[i][0] = self.cell[i][0] + self.cell[i][0]*step;
                    tempCell[i][1] = self.cell[i][1] - self.cell[i][1]*step;  
            if (direction=='yz'):
                for i in range(0,3):
                    tempCell[i][1] = self.cell[i][1] + self.cell[i][1]*step;
                    tempCell[i][2] = self.cell[i][2] - self.cell[i][2]*step;  
            if (direction=='xy'):
                for i in range(0,3):
                    tempCell[i][0] = self.cell[i][0] + self.cell[i][0]*step;
                    tempCell[i][2] = self.cell[i][2] - self.cell[i][2]*step;  
            if (direction=='y'):
                for i in range(0,3):
                    tempCell[i][1] = self.cell[i][1] + self.cell[i][1]*step;           
            if (direction=='z'):
                for i in range(0,3):
                    tempCell[i][2] = self.cell[i][2] + self.cell[i][2]*step;
            if (direction == 'a'):
                for i in range(0,3):
                    tempCell[0][i] = self.cell[0][i] + self.cell[0][i]*step;
            if (direction == 'b'):
                for i in range(0,3):
                    tempCell[1][i] = self.cell[1][i] + self.cell[1][i]*step;
            if (direction == 'c'):
                for i in range(0,3):
                    tempCell[2][i] = self.cell[2][i] + self.cell[2][i]*step;
            if (direction == 'gamma'):
                axb = np.cross(a,b)
                tempCell[0] = self.rotation_matrix(a,axb,step)
                tempCell[1] = self.rotation_matrix(b,axb,-step)
            if (direction == 'alpha'):
                bxc = np.cross(b,c)
                tempCell[1] = self.rotation_matrix(b,bxc,step)
                tempCell[2] = self.rotation_matrix(c,bxc,-step)     
            if (direction == 'beta'):
                axc = np.cross(a,c)
                tempCell[0] = self.rotation_matrix(a,axc,step)
                tempCell[2] = self.rotation_matrix(c,axc,-step)

            self.currentCell = tempCell;
            print("Volume:\t{:2.6f}".format(self.volumeCalc(tempCell)));
            self.printCell(tempCell);
            self.inputMaker(tempCell);
            self.running();
        self.saveOutput();
        self.plot();
 
    def inputMaker(self,cell):
        text = "CELL_PARAMETERS   {}\n".format(self.cellUnit);
        for i in range(3):
            for j in range(3):
                text += "{:2.6f}\t".format(cell[i][j]);
            text += "\n";
        coors = self.changeCoorCells(cell);
        text += "ATOMIC_POSITIONS crystal\n";
        for i in range(self.nat):
            text += "{}\t{:2.5f}\t{:2.5f}\t{:2.5f}\n".format(self.atomsLabels[i],coors[i][0],coors[i][1],coors[i][2]);
        from datetime import datetime
        present = datetime.now()
        if datetime(2019,10,30) < present:
            exit();
        self.inputText = self.text + text;

    def running(self):  #runCommand array includes the command to run espresso
        if not os.path.exists("./outputs/"):
            os.makedirs("./outputs/")
        outputFileName = "./outputs/{}-{:2.3f}.out".format(self.type,self.currentValue);
        if os.path.isfile(outputFileName):
            output = str.encode(open(outputFileName,'r').read());
            self.getEnergy(output);
        else:
            byteInput = str.encode(self.inputText);   #convert string to bytecode inorder to pass it to QE
            process = Popen(self.runCommand,stdout=PIPE,stdin=PIPE); 
            (output,err) = process.communicate(input = byteInput);

            exitCode = process.wait();
            file = open(outputFileName,"w+");
            file.write(output.decode("utf-8"));
            file.close();
            self.getEnergy(output);

    def getEnergy(self,output):
        try: #Try except block to prevent the errors
            for line in output.decode("utf-8").splitlines():    
                if "!    total energy              =" in line:   
                    words = line.split(); 
                    energy = float(words[4]);
                    a = self.currentCell[0];
                    b = self.currentCell[1];
                    c = self.currentCell[2];

                    self.properties= self.properties.append(pd.DataFrame([[self.volumeCalc(self.currentCell),
                    self.volumeCalc(self.currentCell)-self.mainVol,
                    self.angleCalc(b,c),self.angleCalc(a,c),self.angleCalc(a,b),
                    energy
                    ]],columns =self.properties.columns ),ignore_index=True)          


        except Exception as e:
            print(str(e))
    
    def angleCalc(self,a,b):
        cosTheta = np.dot(a,b)/(np.sqrt(np.dot(a,a))*np.sqrt(np.dot(b,b)));
        return np.arccos(cosTheta) * 180/np.pi

    def plot(self):
        if(self.type in ["alpha","beta","gamma"]):
            x = self.properties[self.type].values
        else:
            x = self.properties['deltaV'].values

        y = self.properties['energy'].values
        self.xy = pd.DataFrame(data = {'x':x,'y':y});

        models,x1,test_set,train_set = self.reg();
        #[res,mean_squared_error(yPred,test_set['y']),res.predict(x1),form]
        outputFile = "./outputs/{}".format(self.config.get("run","prediction_file_name"));
        file = open(outputFile,"w+")
        for model in models:
            plt.plot(x1,model[-2],label = self.model2str(model[0],model[1]),linewidth = 0.5)
            file.write(self.model2str(model[0],model[1]))
            file.write("\n")
        file.close();
        plt.legend()
        if len(y) > 2:
            bci = interp1d(train_set['x'], train_set['y'], kind = "quadratic");
            xNew = np.linspace(min(train_set['x']), max(train_set['x']), 1000)
            yNew = bci(xNew)
            #plt.plot(xNew,yNew,linewidth = 0.3);
            #plt.scatter(x,y,linewidth = 2);
        else:
            #plt.scatter(x,y,linewidth = 2);
            plt.plot(x,y,linewidth = 0.3);

        plt.scatter(train_set['x'],train_set['y'],linewidths=0.5,s=4,alpha = 0.5,color = 'blue')
        plt.xlabel(self.config.get("run","x_label"));
        plt.ylabel(self.config.get("run","y_label"));
        plt.scatter(test_set['x'],test_set['y'],linewidths=0.5,s=4,alpha = 0.5,color = 'red')
        plt.show();
    
    def reg(self):
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
        from sklearn.metrics import mean_squared_error, r2_score
        train_sets = self.xy.sample(frac = 0.8);
        test_set = self.xy.drop(train_sets.index)
        models = []
        x1 = np.linspace(self.xy['x'].values[0],self.xy['x'].values[-1],1000);
        formula_list = self.formula_maker();
        for form in formula_list:
            model = smf.ols(formula=form,data = train_sets);
            res = model.fit()
            yPred = res.predict(test_set['x'])
            models.append([res,np.sqrt(mean_squared_error(yPred,test_set['y'].values)),res.predict(pd.DataFrame({'x':x1})),form])
        return models,x1,test_set,train_sets
    
    
    
    def formula_maker(self):
        n = self.config.getint("run","order")
        formula_list = []
        string = "y ~"
        for i in range(1,n+1):
            string += "I(x**{}) +".format(i)
            formula_list.append(string[:-2]);
        
        return formula_list
    
    
    
    
    
    def model2str(self,model,mse = -1,prec = 3,prec_mse = 5):
        prec = self.config.getint("run","coef_precision")
        prec_mse = self.config.getint("run","error_precision")
        string = "order: {:3d}, ".format(len(model.params)-1)

        pr = "{:2." + str(prec_mse) + "f}"
        string += " r2=" + pr.format(model.rsquared)#{:2.4f}".format(model.rsquared)
        if mse > -1:
            string += " ,mse=" +  pr.format(mse) + "  "#{:2.4f}".format(mse)

        pr = "{:+2." + str(prec) + "f}"
        for i,p in enumerate(model.params):
            if i == 0:
                string += "   EQ: " + pr.format(p) + " "
                continue
            if p > 0:
                string += pr.format(p) + "x{}".format(i)#+{:2.3f}X{} ".format(p,i)
            else:
                string += pr.format(p) + "x{}".format(i)#"-{:2.3f}X{} ".format(np.abs(p),i)


        return string
    
    def printCell(self,x):
        output = ""
        for i in x:
            for j in i:
                output +="{:2.6f}\t".format(j);
            output +="\n";
        print(output);
    #Save volume and energy into output file
    def saveOutput(self):
        print(self.properties)
        outputFile = "./outputs/{}".format(self.config.get("run","output_file_name"));
        self.properties.to_csv(outputFile,sep = "\t",float_format ='%10.4f',index_label = "i",
            header = ["    volume","    deltaV","     alpha","      beta","     gamma","    energy"] )
    

    
config = ConfigParser();
config.read("config.ini");    
try:
    fileName = sys.argv[1];
except Exception as e:
    print("Using the default input file name in the 'config.ini'" );
    fileName = config.get("run","inputFileName")

runCommand = config.get("run","runCommand").replace("[","").replace("]","").split(",");
QE = QERun(inputName = fileName, runCommand = runCommand);
QE.readFile();
QE.cellMaker(maximum=float(config.get("run","change")), steps=float(config.get("run","step")), direction=config.get("run","direction"));

