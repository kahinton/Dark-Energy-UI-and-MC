from __future__ import division
from Tkinter import *
import numpy as np
import pylab as py
import os
import sys
from ui_pack import Sing, Show

# Create the class of definitions for the applet
class DE_Mod_Gen(Tk):
    def __init__(self, parent):
        """Initiate the instance of Tkinter"""
        Tk.__init__(self, parent)
        self.parent = parent
        self.initialize()

    def initialize(self):
        """Create the Tkinter objects of DE_Mod_Gen"""
        self.grid()

        #Create the Home Screen of the applet with buttons to access the other pages        
        
        self.Home = Frame(self, height = 475, width = 700)
        self.Home.grid(column = 0, row = 0, columnspan = 3, sticky = "N")
        self.Home.grid_propagate(0)
        
        self.IntroStr = StringVar()
        self.Welcome = Label(self.Home, textvariable = self.IntroStr, anchor = "center", fg = "black")
        self.Welcome.grid(column = 0, row = 0, columnspan = 4)
        self.IntroStr.set(u"Welcome to DE_Mod_Gen! Using the buttons below you may choose to run a single instance of any")
        
        self.IntroStr2 = StringVar()
        self.Welcome2 = Label(self.Home, textvariable = self.IntroStr2, anchor = "center", fg = "black")
        self.Welcome2.grid(column = 0, row = 1, columnspan = 4, pady = 3)
        self.IntroStr2.set(u"model previously run with the Monte Carlo generator, or analyze the bulk results from any of those")
        
        self.IntroStr3 = StringVar()
        self.Welcome3 = Label(self.Home, textvariable = self.IntroStr3, anchor = "center", fg = "black")
        self.Welcome3.grid(column = 0, row = 2, columnspan = 4, pady = 3)
        self.IntroStr3.set(u"tests. You can also open the Monte Carlo code if you have the appropriate program.")
        
        self.SinglePage = Button(self.Home, text = u"Run a Single Model", command = self.Single1, width = 35)
        self.SinglePage.grid(column = 2, row = 3, padx = 10, pady = 30)
        
        self.ResultPage = Button(self.Home, text = u"Results from Monte Carlo Generator", command = self.Result1, width = 35)
        self.ResultPage.grid(column = 2, row = 4, padx = 10, pady = 30)
                
        self.MontePage = Button(self.Home, text = u"Open Monte Carlo Code (requires iPython)", command = self.MonteOpen, width = 35)
        self.MontePage.grid(column = 2, row = 5, padx = 10, pady = 30)        
        
        # Create window for running a single model
      
        self.SingleInput = Frame(self, height = 475, width = 700)
        self.SingleInput.grid(column = 0, row = 0, columnspan = 3, sticky = "N")
        self.SingleInput.grid_propagate(0)
        
        # Find All Usable Models
        
        __location__ = os.path.dirname(sys.argv[0]) 
        loc = os.path.join(__location__, 'Parameters/')
        self.Model_Files = os.listdir(loc)
        self.Model_Names = []
        for j in self.Model_Files:
            q = os.path.splitext(j)
            name = str(q[0])
            self.Model_Names.append(name)
        self.SingleModsVar = StringVar()
        self.StringModsOpt = OptionMenu(self.SingleInput, self.SingleModsVar, *self.Model_Names)
        self.StringModsOpt.grid(column = 2, row = 0, pady = 5)
        self.SingleModsVar.set('Click')
        
        # Create All Labels in the single model window
        
        self.Labels_Name = ['Mod', 'Input', 'Value', 'Suggested','x1Disc','x1Ran','x2Disc','x2Ran', 'y1Disc','y1Ran','y2Disc','y2Ran','c1Disc','c1Ran','c2Disc','c2Ran','uDisc','uRan','nDisc','nRan','xaxis','yaxis']
        self.Labels_Grid = np.array([[0,0,2],[0,2,2],[2,2,1],[3,2,1],[0,4,2],[3,4,1],[0,5,2],[3,5,1],[0,6,2],[3,6,1],[0,7,2],[3,7,1],[0,8,2],[3,8,1],[0,9,2],[3,9,1],[0,10,2],[3,10,1],[0,11,2],[3,11,2],[4,2,2],[4,6,2]])
        self.Labels_Vars = ['Please Select a Model ->','Model Parameters','Enter Values','Suggested',"x1 (Early Field Kinetic Energy)","(---)","x2 (Late Field Kinetic Energy)","(---)","y1 (Early Field Potential Energy)","(---)","y2 (Late Field Potential Energy)","(---)","c1 (Early Field Potential Coefficient)","(---)","c2 (Late Field Potential Coefficient)","(---)","u (Describes Radiation Density)","(---)","Total Number of Fields","---",'Pick an x axis','Pick Data to Plot']
                    
        for z in xrange(len(self.Labels_Name)):
            exec('self.Single'+str(self.Labels_Name[z])+'Var = StringVar()')
            exec('self.Single'+str(self.Labels_Name[z])+'Lab = Label(self.SingleInput, textvariable = self.Single'+str(self.Labels_Name[z])+'Var, anchor = "w", fg = "black")')
            exec('self.Single'+str(self.Labels_Name[z])+'Lab.grid(column = '+str(self.Labels_Grid[z,0])+',row = '+str(self.Labels_Grid[z,1])+', columnspan = '+str(self.Labels_Grid[z,2])+', sticky = "EW", pady = 5)')
            exec('self.Single'+str(self.Labels_Name[z])+'Var.set(u"'+str(self.Labels_Vars[z])+'")')
        
        # Create all checkboxes in the single model window
        
        self.Box_Name = ['a','z','OM','ODE','OR','W','D']
        self.Box_Grid = np.array([[4,4],[4,5],[4,7],[4,8],[4,9],[4,10],[4,11]])                   
        self.Box_Text = ['a','z','Omega_M','Omega_DE','Omega_Rad','w_DE','Growth D/a']
        
        for q in xrange(len(self.Box_Name)):
            exec('self.Single'+str(self.Box_Name[q])+'Var = StringVar()')
            exec('self.Single'+str(self.Box_Name[q])+'Box = Checkbutton(self.SingleInput, variable =self.Single'+str(self.Box_Name[q])+'Var, text = "'+str(self.Box_Text[q])+'")')
            exec('self.Single'+str(self.Box_Name[q])+'Box.grid(column = '+str(self.Box_Grid[q,0])+', row = '+str(self.Box_Grid[q,1])+', sticky = "W", pady = 2)')
            exec('self.Single'+str(self.Box_Name[q])+'Var.set(0)')
        
        # Create all entry objects in the single model window
        
        self.Entry_Name = ['x1','x2','y1','y2','c1','c2','u','n']
        self.Entry_Grid = np.array([[2,4],[2,5],[2,6],[2,7],[2,8],[2,9],[2,10],[2,11]])
                                    
        for p in xrange(len(self.Entry_Name)):
            exec('self.Single'+str(self.Entry_Name[p])+'Var = StringVar()')
            exec('self.Single'+str(self.Entry_Name[p])+'Ent = Entry(self.SingleInput, textvariable = self.Single'+str(self.Entry_Name[p])+'Var)')
            exec('self.Single'+str(self.Entry_Name[p])+'Ent.grid(column = '+str(self.Entry_Grid[p,0])+', row = '+str(self.Entry_Grid[p,1])+', sticky = "EW", pady = 5)')
            exec('self.Single'+str(self.Entry_Name[p])+'Var.set(u"---")')
            
        # Buttons to initiate the model with users parameters or to go back to the Home window
        
        self.SingleSelect = Button(self.SingleInput, text = u'Then Click Here', command = self.SingleCheck)
        self.SingleSelect.grid(column = 3, row = 0, pady = 5)
        
        self.SingleStart = Button(self.SingleInput, text = u"Run Model", command = self.SingleRun)
        self.SingleStart.grid(column = 2, row = 12 )
        
        self.SingleGoHome = Button(self.SingleInput, text = u"Home", command = self.BackToHome)
        self.SingleGoHome.grid(column = 3, row = 12)
        
        self.SingleInput.grid_remove()

        # Create the page to access and view results of the Monte Carlo generator       
        
        self.Result = Frame(self, height = 475, width = 700)
        self.Result.grid(column = 0, row = 0, sticky = "N")
        self.Result.grid_propagate(0)
        
        # Create the labels in the Result window
        
        self.RLabel_Name = ['Exp1', 'Exp2','Exp3','Mod1','X1','Y1','Mod2','X2','Y2']
        self.RLabel_Grid = np.array([[0,0,5],[0,1,5],[0,2,5],[0,3,1],[0,5,1],[0,6,1],[3,3,1],[3,5,1],[3,6,1]])
        self.RLabel_Vars = ['From this window you can analyze the results from the Monte Carlo tests of your choice.','The selected parameters can be plotted in (x1,y1) and (x2,y2) planes. Available parameters are', 'w_0, w_a, w_p, w_DE, Omega_M, Omega_DE, x1, x2, y1, y2, c1, c2, u, Lambda_1, Lambda_2, and n','First Model','Graph 1 x-axis','Graph 1 y-axis','Second Model','Graph 2 x-axis','Graph 2 y-axis']
        
        for b in xrange(len(self.RLabel_Name)):
            exec('self.Result'+str(self.RLabel_Name[b])+'Var = StringVar()')
            exec('self.Result'+str(self.RLabel_Name[b])+'Lab = Label(self.Result, textvariable = self.Result'+str(self.RLabel_Name[b])+'Var, anchor = "center")')
            exec('self.Result'+str(self.RLabel_Name[b])+'Lab.grid(column = '+str(self.RLabel_Grid[b,0])+',row = '+str(self.RLabel_Grid[b,1])+', columnspan = '+str(self.RLabel_Grid[b,2])+', sticky = "EW", pady = 5)')
            exec('self.Result'+str(self.RLabel_Name[b])+'Var.set(u"'+str(self.RLabel_Vars[b])+'")')
        
        # Create the Optionmenus in the Result window
        
        self.ROption_Name = ['Pick1','Pick2']
        self.ROption_Grid = np.array([[1,3],[4,3]])
        Rloc = os.path.join(__location__, 'Data/')
        self.RModel_Files = os.listdir(Rloc)
        self.ROption_List = []
        for j in self.RModel_Files:
            Rq = os.path.splitext(j)
            Rname = str(Rq[0])
            self.ROption_List.append([Rname])
        
        for r in xrange(len(self.ROption_Name)):
            exec('self.Result'+str(self.ROption_Name[r])+'Var = StringVar()')
            exec('self.Result'+str(self.ROption_Name[r])+'Opt = OptionMenu(self.Result, self.Result'+str(self.ROption_Name[r])+'Var,*'+str(self.ROption_List)+')')
            exec('self.Result'+str(self.ROption_Name[r])+'Opt.config(width = 8)')
            exec('self.Result'+str(self.ROption_Name[r])+'Opt.grid(column = '+str(self.ROption_Grid[r,0])+', row = '+str(self.ROption_Grid[r,1])+', pady = 10)')
            exec('self.Result'+str(self.ROption_Name[r])+'Var.set(u"Click")')
         
         # Create the Entry objects in the Result window   
          
        self.REntry_Name = ['XAxis1', 'YAxis1','XAxis2','YAxis2']
        self.REntry_Exam = ['e.g. w_0', 'e.g. w_a', 'e.g. np.log10(x1)', 'e.g. x1/y1']
        self.REntry_Grid = np.array([[1,5],[1,6],[4,5],[4,6]])
        
        for k in xrange(len(self.REntry_Name)):
            exec('self.Result'+str(self.REntry_Name[k])+'Var = StringVar()')
            exec('self.Result'+str(self.REntry_Name[k])+'Ent = Entry(self.Result, textvariable = self.Result'+str(self.REntry_Name[k])+'Var)')
            exec('self.Result'+str(self.REntry_Name[k])+'Ent.grid(column = '+str(self.REntry_Grid[k,0])+', row = '+str(self.REntry_Grid[k,1])+', sticky = "EW", pady = 10)')
            exec('self.Result'+str(self.REntry_Name[k])+'Var.set(u"'+self.REntry_Exam[k]+'")')
        
        # Bottons to get the results of the query or to go back to the Home window
        
        self.ResultGet = Button(self.Result, text = u"Graph Results", command = self.ResultView)
        self.ResultGet.grid(column = 2, row = 7, pady = 5)
        
        self.ResultGoHome = Button(self.Result, text = u"Home", command = self.BackToHome)
        self.ResultGoHome.grid(column = 2, row = 8, pady = 5)
        
        self.Result.grid_remove()

        # Create a status bar to inform the user of the applets actions
        
        self.statusVariable = StringVar()
        self.status = Label(self, textvariable = self.statusVariable, anchor = "w", fg = "yellow", bg = "blue")
        self.status.grid(column = 0, row = 10,columnspan = 5, sticky = 'EWS')
        self.statusVariable.set(u"Welcome to DE_Mod_Gen")

        # Configure how the applet can be resized by the user
        
        self.grid_columnconfigure(0, weight = 1)
        self.resizable(True, False)
        self.update()
        self.geometry(self.geometry())       
         
    def Single1(self):
        """Single1 takes the user to the Single Model window"""
        self.statusVariable.set("Run a single model with your chosen parameters")
        self.Home.grid_remove()
        self.Result.grid_remove()
        self.SingleInput.grid()
        self.Singlex1Ent.focus_set()
        self.Singlex1Ent.selection_range(0, END)
        
    def SingleCheck(self):
        """Set the single model inputs to fiducial values for Quintessence"""
        
        __location__ = os.path.dirname(sys.argv[0]) 
        loc = os.path.join(__location__, 'Parameters/')
        SingleSetList = ['x1','x2','y1','y2','c1','c2','u','n','x1Ran','x2Ran','y1Ran','y2Ran','c1Ran','c2Ran','uRan','nRan']
        ModSet = np.genfromtxt(loc+str(self.SingleModsVar.get())+'.dat', delimiter = 'newcolumn')
        SingleSetVars = [str((ModSet[5]+ModSet[6])/2.0),str((ModSet[7]+ModSet[8])/2.0),str((ModSet[9]+ModSet[10])/2.0),str((ModSet[11]+ModSet[12])/2.0),str((ModSet[13]+ModSet[14])/2.0),str((ModSet[15]+ModSet[16])/2.0),str((ModSet[17]+ModSet[18])/2.0),int((ModSet[19])),
                        "u'("+str(ModSet[5])+","+str(ModSet[6])+")'","u'("+str(ModSet[7])+","+str(ModSet[8])+")'","u'("+str(ModSet[9])+","+str(ModSet[10])+")'","u'("+str(ModSet[11])+","+str(ModSet[12])+")'",
                        "u'("+str(ModSet[13])+","+str(ModSet[14])+")'","u'("+str(ModSet[15])+","+str(ModSet[16])+")'","u'("+str(ModSet[17])+","+str(ModSet[18])+")'","u'Integer'"]
        
        for v in xrange(len(SingleSetList)):
            exec('self.Single'+str(SingleSetList[v])+'Var.set('+str(SingleSetVars[v])+')')
        
        self.Singlex1Ent.focus_set()
        self.Singlex1Ent.selection_range(0, END)

     
    def MonteOpen(self):
        """Opens the file to run the Monte Carlo generator which provides information to the UI"""
        __location__ = os.path.dirname(sys.argv[0]) 
        os.system(__location__+'/Monte.py') 
        self.statusVariable.set("Please restart UI after running the Monte Carlo generator")  
         
    def Result1(self):
        """Result1 takes the user to the Result window"""
        self.statusVariable.set("Results from Monte Carlo Generator")
        self.Home.grid_remove()
        self.SingleInput.grid_remove()
        self.Result.grid()
        self.ResultXAxis1Ent.focus_set()
        self.ResultXAxis1Ent.selection_range(0, END)
      
    def BackToHome(self):
        """BackToHome takes the user back to the Home Window"""
        self.SingleInput.grid_remove()
        self.Result.grid_remove()
        self.Home.grid()
        self.statusVariable.set("Welcome to DE_Mod_Gen")
        
    def GraphMaker(self,name,c,d,e):
        """Produces graphs based on the user entered equations in the Result Window"""
        
        model = name.rstrip("']")
        model = model.lstrip("'[")
        model = str(model)
        name = ''
        for i in xrange(len(model)):
            if model[i] == '_':
                name += ' '
            else:
                name += model[i]
                
        __location__ = os.path.dirname(sys.argv[0]) 
        loc = os.path.join(__location__, 'Data/')
        
        Params = ["w_0", "w_a", "w_p", "w_DE", "Omega_M", "Omega_DE", "x1", "x2", "y1", "y2", "c1", "c2", "u", "Lambda_1", "Lambda_2", "n"]
        for entry in xrange(len(Params)):
            exec( Params[entry] + " = np.genfromtxt(loc+str(model)+'.dat', usecols = "+str(entry)+", skip_header = 1)" )
                
        A = eval(c)
        B = eval(d)
            
        py.figure(int(e))
        py.plot(A, B, 'b.')
        py.xlabel(c, fontsize = 55)
        py.ylabel(d, fontsize = 55)
        py.title('Results from '+str(name), fontsize = 55)
        py.tick_params(labelsize = 35, size = 15, width = 5, top = 0, right = 0)
            
        py.show()
        
    def ResultView(self):
        """ResultView shows the results of Monte Carlo tests for the given user selected models"""
        
        self.statusVariable.set("Gathering Data")
        if str(self.ResultPick1Var.get()) != 'Click':
            if str(self.ResultXAxis1Var.get()) and str(self.ResultXAxis1Var.get()) != 'Python Input':
                self.GraphMaker(self.ResultPick1Var.get(),self.ResultXAxis1Var.get(),self.ResultYAxis1Var.get(),1)
        if str(self.ResultPick2Var.get()) != 'Click':
            if str(self.ResultXAxis2Var.get()) and str(self.ResultXAxis2Var.get()) != 'Python Input':
                self.GraphMaker(self.ResultPick2Var.get(),self.ResultXAxis2Var.get(),self.ResultYAxis2Var.get(),2)       
        self.statusVariable.set("Displaying Data")
 
    def SingleRun(self):
        """SingleRun takes users input parameters and runs a single instance of the equations of evolution"""
        
        self.statusVariable.set( "-Running Simulation-" )
        __location__ = os.path.dirname(sys.argv[0]) 
        loc = os.path.join(__location__, 'Parameters/')
        
        Model = np.genfromtxt(loc+str(self.SingleModsVar.get())+'.dat', delimiter = 'newcolumn', dtype = 'string')
        x1 = eval(self.Singlex1Ent.get())
        x2 = eval(self.Singlex2Ent.get())       
        y1 = eval(self.Singley1Ent.get())
        y2 = (eval(self.Singley2Ent.get()))
        u = eval(self.SingleuEnt.get())
        c1 = eval(self.Singlec1Ent.get())
        c2 = eval(self.Singlec2Ent.get())
        n = int(self.SinglenEnt.get())-1
        
        A,Z,OMat,ODE,ORad,WPHI,D = Sing(x1,x2,y1,y2,c1,c2,u,n,Model)
        
        axis1 = int(self.SingleaVar.get())
        axis2 = int(self.SinglezVar.get())
        Check1 = self.SingleOMVar.get()
        Check2 = self.SingleODEVar.get()
        Check3 = self.SingleORVar.get()
        Check4 = self.SingleWVar.get()
        Check5 = self.SingleDVar.get()
        
        Out = Show(axis1, axis2, Check1, Check2, Check3, Check4, Check5, A, Z, OMat, ODE, ORad, WPHI, D, Model, n)
        
        self.statusVariable.set(Out)
        
# Run an instance of the applet

if __name__ == "__main__":
    app = DE_Mod_Gen(None)
    app.title('DE_Mod_Gen')
    app.mainloop()