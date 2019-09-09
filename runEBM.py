from subprocess import call
import os


os.chdir("./model")
call("./driver",shell=True) #run the driver
    

