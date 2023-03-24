import sys
import os
import readline

readline.parse_and_bind('tab: complete')
readline.parse_and_bind('set editing-mode vi')

binpath='/home/nico/Workspace/GeneralFanoCICode/PythonScript/'

print "" 
print "This is a Python script that will help you to prepare the input files you need to"
print "1) Run CI calculations"
print "2) Compute total decay widths"
print "3) Compute partial decay widths"
print "" 

opt = int(raw_input('What do you want to do?  '))

if opt==1:         # CI calculations
	print "" 
	nxml = int(raw_input('How many XML files do you have?  '))
	fxml = []
	sedopt = []
	for i in range(nxml):
		print ""
		iopt = raw_input('Is it for a single CI calc.? y/Y/yes/YES or n/N/no/NO  ')
		if(iopt=='y' or iopt=='Y' or iopt=='yes' or iopt=='YES'):
			fxml.append( raw_input('Enter the XML filename  ') )
			sedopt.append([-1,-1])
		else:
			fxml.append( raw_input('Enter the XML filename  ') )
			imin = int(raw_input('Enter the index of the first orbital  '))
			imax = int(raw_input('Enter the index of the last orbital  '))
			sedopt.append([imin,imax])

	print ""
	print "A shell script will be generated"
        scriptname=raw_input('Enter its filename ')
	fshell = open(scriptname,"w")
	ninput = -1
	for i in range(nxml):
	    if(sedopt[i][0]==-1):
		ninput+=1
		print >> fshell, "python ", binpath+"ormas_like.py ", fxml[i], "> a"
 		print >> fshell, "cat header.txt a > ", "input"+str(ninput)
	    else:
		for iorb in range(sedopt[i][0],sedopt[i][1]):
			ninput+=1
			print >> fshell, 'sed "s/X/'+str(iorb)+'/g"', fxml[i], "> a"
			print >> fshell, "python ", binpath+"ormas_like.py a > b"
 			print >> fshell, "cat header.txt b > ", "input"+str(ninput)
	
	print "input_myci will be generated"
	fci =  open("input_myci","w")
	print >> fci, ninput+1
	for i in range(ninput+1):
		print >> fci, 'input'+str(i)

        cmd="chmod u+x "+scriptname
        os.system(cmd)
	print >> fshell, "rm a b header.txt"
	print "Dont forget to run ./"+scriptname
	print "Finally Run myCI < input_myci"
elif opt==2:         # Total widths calculations
	print "" 
	fcidecay = []
	ncista = []
	ncidecay = int(raw_input('How many CI files do you have for the decaying states?  '))
	if(ncidecay==1):
		fcidecay.append( raw_input('Enter the CI filename  ') )
		ncista.append(int(raw_input('How many CI states do you have?  ') ))
	else:	
		dopt = int(raw_input('Do you want to enter them manually (1) or through a loop (2)?'))
		if (dopt==1):
			for i in range(ncidecay):
				fcidecay.append( raw_input('Enter the CI filename  ') )
				ncista.append(int(raw_input('How many CI states do you have?  ') ))
		else:
			print "The files should be named like ftemplate_index"
			print "and you should have the same number of CI states in each file"
			ftemp = raw_input('Enter the filename template  ')
			imin = int(raw_input('Enter the first index  '))
			ista = int(raw_input('How many CI states do you have?  '))
			for i in range(ncidecay):
				fcidecay.append(ftemp+str(imin+i))
				ncista.append(ista)

	print "" 
	fcifinal = []
	ncifsta = []
	sedopt = []
	ncifinal = int(raw_input('How many CI files do you have for the final states?  '))
	if(ncifinal==1):
		fcifinal.append( raw_input('Enter the CI filename  ') )
		ncifsta.append(int(raw_input('How many CI states do you have?  ') ))
	else:	
		dopt = int(raw_input('Do you want to enter them manually (1) or through a loop (2)?'))
		if (dopt==1):
			for i in range(ncifinal):
				fcifinal.append( raw_input('Enter the CI filename  ') )
				ncifsta.append(int(raw_input('How many CI states do you have?  ') ))
		else:
			print "The files should be named like ftemplate_index"
			print "and you should have the same number of CI states in each file"
			ftemp = raw_input('Enter the filename template  ')
			imin = int(raw_input('Enter the first index  '))
			ista = int(raw_input('How many CI states do you have?  '))
			for i in range(ncifinal):
				fcifinal.append(ftemp+str(imin+i))
				ncifsta.append(ista)

	fname=raw_input('Enter the Fano input filename that will be generated ')
	fanoinp=open(fname,"w")
		
	k=-1
	print >> fanoinp, ncidecay*ncifinal
	for i in range(ncidecay):
		for j in range(ncifinal):
			k+=1
			print >> fanoinp, ncista[i],ncifsta[j],fcidecay[i],fcifinal[j],' fano'+str(k)	

	print "" 
	print "" 
	print "Run Fano < ",fname
elif opt==3:         # Partial widths calculations
	print "" 
	fcidecay = []
	ncista = []
	ncidecay = int(raw_input('How many CI files do you have for the decaying states?  '))
	if(ncidecay==1):
		fcidecay.append( raw_input('Enter the CI filename  ') )
		ncista.append(int(raw_input('How many CI states do you have?  ') ))
	else:	
		dopt = int(raw_input('Do you want to enter them manually (1) or through a loop (2)?'))
		if (dopt==1):
			for i in range(ncidecay):
				fcidecay.append( raw_input('Enter the CI filename  ') )
				ncista.append(int(raw_input('How many CI states do you have?  ') ))
		else:
			print "The files should be named like ftemplate_index"
			print "and you should have the same number of CI states in each file"
			ftemp = raw_input('Enter the filename template  ')
			imin = int(raw_input('Enter the first index  '))
			ista = int(raw_input('How many CI states do you have?  '))
			for i in range(ncidecay):
				fcidecay.append(ftemp+str(imin+i))
				ncista.append(ista)

        print ""
        fcitarg = []
        ncitsta = []
        ncitarg = int(raw_input('How many CI files do you have for the target states?  '))
        if(ncitarg==1):
                fcitarg.append( raw_input('Enter the CI filename  ') )
                ncitsta.append(int(raw_input('How many CI states do you have?  ') ))
        else:
                dopt = int(raw_input('Do you want to enter them manually (1) or through a loop (2)?'))
                if (dopt==1):
                        for i in range(ncitarg):
                                fcitarg.append( raw_input('Enter the CI filename  ') )
                                ncitsta.append(int(raw_input('How many CI states do you have?  ') ))
                else:
                        print "The files should be named like ftemplate_index"
                        print "and you should have the same number of CI states in each file"
                        ftemp = raw_input('Enter the filename template  ')
                        imin = int(raw_input('Enter the first index  '))
                        ista = int(raw_input('How many CI states do you have?  '))
                        for i in range(ncitarg):
                                fcitarg.append(ftemp+str(imin+i))
                                ncitsta.append(ista)


	print "" 
	fcifinal = []
	ncifsta = []
	sedopt = []
	ncifinal = int(raw_input('How many CI files do you have for the final states?  '))
	if(ncifinal==1):
		fcifinal.append( raw_input('Enter the CI filename  ') )
		ncifsta.append(int(raw_input('How many CI states do you have?  ') ))
	else:	
		dopt = int(raw_input('Do you want to enter them manually (1) or through a loop (2)?'))
		if (dopt==1):
			for i in range(ncifinal):
				fcifinal.append( raw_input('Enter the CI filename  ') )
				ncifsta.append(int(raw_input('How many CI states do you have?  ') ))
		else:
			print "The files should be named like ftemplate_index"
			print "and you should have the same number of CI states in each file"
			ftemp = raw_input('Enter the filename template  ')
			imin = int(raw_input('Enter the first index  '))
			ista = int(raw_input('How many CI states do you have?  '))
			for i in range(ncifinal):
				fcifinal.append(ftemp+str(imin+i))
				ncifsta.append(ista)

	fname=raw_input('Enter the Fano input filename that will be generated ')
	fanoinp=open(fname,"w")
		
	k=-1
	print >> fanoinp, ncidecay*ncifinal*ncitarg
	for i in range(ncidecay):
	    for it in range(ncitarg):
		for j in range(ncifinal):
			k+=1
			print >> fanoinp, ncista[i],ncifsta[j],ncitsta[it],fcidecay[i],fcifinal[j],fcitarg[it],' fanopart'+str(k)	

	print "" 
	print "" 
	print "Run Fano < ",fname
	

	

