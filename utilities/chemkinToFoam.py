import re

# Reactions
numElements = 0
numSpecies = 0
elements = []
species = []
with open('chem.inp') as f:
    for line in f:
        if 'ELEMENTS' in line:
            break
    for line in f:
        line = line.strip()
        if 'END' in line:
            break
        else:
            line = re.split(r'\s+',line)
            for name in line:
                elements.append(name)
                numElements += 1
    
    for line in f:
        if 'SPECIE' in line:
            break
    for line in f:
        line = line.strip()
        if 'END' in line:
            break
        else:
            line = re.split(r'\s+', line)
            for name in line:
                species.append(name)
                numSpecies += 1

with open('reactions','w+') as f:
    f.write('elements\n')
    f.write(f'{numElements}')
    f.write('\n(\n')
    for elem in elements:
        f.write(elem)
        f.write('\n')

    f.write(')\n;\n\nspecies\n')
    f.write(f'{numSpecies}')
    f.write('\n(\n')
    for spec in species:
        f.write(spec)
        f.write('\n')

    f.write(')\n;\n\n')
    f.write('reactions\n{')

EToTa = 4184.0/8314.0 # E:[cal/mol] R:[J/kmol K]
orderDict = {0:1e3, 1:1.0, 2:1e-3, 3:1e-6, 4:1e-9}
numReactions = 0
numbers = '0123456789'
def thirdBodyArrheniusParse(line,f,parser,il):
    f.write('\t\treaction\t\t\"')
    reactionOrder = 0
    reactSpe = line[0].split(parser)
    lhs = reactSpe[0]
    rhs = reactSpe[1]
    reactants = lhs.split('+')
    products = rhs.split('+')
    for j,ireactant in enumerate(reactants):
        if ireactant[0] not in numbers:
            reactionOrder+=1
        else:
            reactionOrder+=int(ireactant[0])
        if j!= (len(reactants)-1):
            f.write(ireactant+' ')
        else:
            pass
        if j < (len(reactants)-2):
            f.write('+ ')
        else:
            pass
    f.write('=')
    for j,iproduct in enumerate(products):
        if j!= (len(products)-1):
            f.write(' '+iproduct)
        else:
            pass
        if j < (len(products)-2):
            f.write(' +')
        else:
            pass
    f.write('";\n')
    f.write('\t\tA\t\t\t\t' + f'{orderDict[reactionOrder]*float(line[1])}' + ';\n')
    f.write('\t\tbeta\t\t\t' + line[2] + ';\n')
    f.write('\t\tTa\t\t\t\t' + f'{float(line[3])*EToTa}' + ';\n')
    f.write('\t\tcoeffs\n' + f'{numSpecies}' + '\n(\n')
    line2 = lines[il+1]
    if '/' in line2:
        line2 = line2.strip()
        line2 = re.split(r'\s+', line2)
        for iLine2 in range(len(line2)):
            line2[iLine2] = line2[iLine2].strip('/')
        speciesCoeffs = []
        for ispecie in range(numSpecies):
            speciesCoeffs.append('1')
            for iLine2,modSpec in enumerate(line2):
                if modSpec == species[ispecie]:
                    speciesCoeffs[ispecie] = line2[iLine2+1]
                    break
                else:
                    continue
    else:
        speciesCoeffs = []
        for ispecie in range(numSpecies):
            speciesCoeffs.append('1')
    for ispecie in range(numSpecies):
        f.write('(' + species[ispecie] + ' ' + speciesCoeffs[ispecie] + ')\n')
    f.write(')\n;\n\t}')

def troeReactionParse(line,f,parser,il):
    f.write('\t\treaction\t\t\"')
    reactionOrder = 0
    reactSpe = line[0].split(parser)
    lhs = reactSpe[0]
    rhs = reactSpe[1]
    lhs = lhs.strip('(+M)')
    rhs = rhs.strip('(+M)')
    reactants = lhs.split('+')
    products = rhs.split('+')
    for j,ireactant in enumerate(reactants):
        if ireactant[0] not in numbers:
            reactionOrder+=1
        else:
            reactionOrder+=int(ireactant[0])
        f.write(ireactant+' ')
        if j != (len(reactants)-1):
            f.write('+ ')
        else:
            pass
    f.write('=')
    for j,iproduct in enumerate(products):
        f.write(' '+iproduct)
        if j != (len(products)-1):
            f.write(' +')
        else:
            pass
    lineLow = lines[il+1]
    lineLow = lineLow.strip()
    lineLow = re.split(r'\s+', lineLow)
    for iLineLow in range(len(lineLow)):
        lineLow[iLineLow] = lineLow[iLineLow].strip('/')
    lineTroe = lines[il+2]
    lineTroe = lineTroe.strip()
    lineTroe = re.split(r'\s+', lineTroe)
    for iLineTroe in range(len(lineTroe)):
        lineTroe[iLineTroe] = lineTroe[iLineTroe].strip('/')
    lineTroe.append('1e-30')
    f.write('";\n')
    f.write('\t\tk0\n\t\t{\n')
    f.write('\t\t\tA\t\t\t' + f'{orderDict[reactionOrder+1]*float(lineLow[1])}' + ';\n')
    f.write('\t\t\tbeta\t\t' + lineLow[2] + ';\n')
    f.write('\t\t\tTa\t\t\t' + f'{float(lineLow[3])*EToTa}' + ';\n\t\t}\n')
    f.write('\t\tkInf\n\t\t{\n')
    f.write('\t\t\tA\t\t\t' + f'{orderDict[reactionOrder]*float(line[1])}' + ';\n')
    f.write('\t\t\tbeta\t\t' + line[2] + ';\n')
    f.write('\t\t\tTa\t\t\t' + f'{float(line[3])*EToTa}' + ';\n\t\t}\n')
    
    f.write('\t\tF\n\t\t{\n')
    f.write('\t\t\talpha\t\t' + lineTroe[1] + ';\n')
    f.write('\t\t\tTsss\t\t' + lineTroe[2] + ';\n')
    f.write('\t\t\tTs\t\t\t' + lineTroe[3] + ';\n')
    f.write('\t\t\tTss\t\t\t' + lineTroe[4] + ';\n\t\t}\n')

    f.write('\t\tthirdBodyEfficiencies\n\t\t{\n')
    f.write('\t\t\tcoeffs\n' + f'{numSpecies}' + '\n(\n')
    line2 = lines[il+3]
    if '/' in line2:
        line2 = line2.strip()
        line2 = re.split(r'\s+', line2)
        for iLine2 in range(len(line2)):
            line2[iLine2] = line2[iLine2].strip('/')
        speciesCoeffs = []
        for ispecie in range(numSpecies):
            speciesCoeffs.append('1')
            for iLine2,modSpec in enumerate(line2):
                if modSpec == species[ispecie]:
                    speciesCoeffs[ispecie] = line2[iLine2+1]
                    break
                else:
                    continue
    else:
        speciesCoeffs = []
        for ispecie in range(numSpecies):
            speciesCoeffs.append('1')
    for ispecie in range(numSpecies):
        f.write('(' + species[ispecie] + ' ' + speciesCoeffs[ispecie] + ')\n')
    f.write(')\n;\n\t\t}\n\t}')

def sriReactionParse(line,f,parser,il):
    f.write('\t\treaction\t\t\"')
    reactionOrder = 0
    reactSpe = line[0].split(parser)
    lhs = reactSpe[0]
    rhs = reactSpe[1]
    lhs = lhs.strip('(+M)')
    rhs = rhs.strip('(+M)')
    reactants = lhs.split('+')
    products = rhs.split('+')
    for j,ireactant in enumerate(reactants):
        if ireactant[0] not in numbers:
            reactionOrder+=1
        else:
            reactionOrder+=int(ireactant[0])
        f.write(ireactant+' ')
        if j != (len(reactants)-1):
            f.write('+ ')
        else:
            pass
    f.write('=')
    for j,iproduct in enumerate(products):
        f.write(' '+iproduct)
        if j != (len(products)-1):
            f.write(' +')
        else:
            pass
    lineLow = lines[il+1]
    lineLow = lineLow.strip()
    lineLow = re.split(r'\s+', lineLow)
    for iLineLow in range(len(lineLow)):
        lineLow[iLineLow] = lineLow[iLineLow].strip('/')
    lineSRI = lines[il+2]
    lineSRI = lineSRI.strip()
    lineSRI = re.split(r'\s+', lineSRI)
    for ilineSRI in range(len(lineSRI)):
        lineSRI[ilineSRI] = lineSRI[ilineSRI].strip('/')
    f.write('";\n')
    f.write('\t\tk0\n\t\t{\n')
    f.write('\t\t\tA\t\t\t' + f'{orderDict[reactionOrder+1]*float(lineLow[1])}' + ';\n')
    f.write('\t\t\tbeta\t\t' + lineLow[2] + ';\n')
    f.write('\t\t\tTa\t\t\t' + f'{float(lineLow[3])*EToTa}' + ';\n\t\t}\n')
    f.write('\t\tkInf\n\t\t{\n')
    f.write('\t\t\tA\t\t\t' + f'{orderDict[reactionOrder]*float(line[1])}' + ';\n')
    f.write('\t\t\tbeta\t\t' + line[2] + ';\n')
    f.write('\t\t\tTa\t\t\t' + f'{float(line[3])*EToTa}' + ';\n\t\t}\n')
    
    f.write('\t\tF\n\t\t{\n')
    f.write('\t\t\ta\t\t\t' + lineSRI[1] + ';\n')
    f.write('\t\t\tb\t\t\t' + lineSRI[2] + ';\n')
    f.write('\t\t\tc\t\t\t' + lineSRI[3] + ';\n')
    f.write('\t\t\td\t\t\t' + lineSRI[4] + ';\n')
    f.write('\t\t\te\t\t\t' + lineSRI[5] + ';\n\t\t}\n')

    f.write('\t\tthirdBodyEfficiencies\n\t\t{\n')
    f.write('\t\t\tcoeffs\n' + f'{numSpecies}' + '\n(\n')
    line2 = lines[il+3]
    if '/' in line2:
        line2 = line2.strip()
        line2 = re.split(r'\s+', line2)
        for iLine2 in range(len(line2)):
            line2[iLine2] = line2[iLine2].strip('/')
        speciesCoeffs = []
        for ispecie in range(numSpecies):
            speciesCoeffs.append('1')
            for iLine2,modSpec in enumerate(line2):
                if modSpec == species[ispecie]:
                    speciesCoeffs[ispecie] = line2[iLine2+1]
                    break
                else:
                    continue
    else:
        speciesCoeffs = []
        for ispecie in range(numSpecies):
            speciesCoeffs.append('1')
    for ispecie in range(numSpecies):
        f.write('(' + species[ispecie] + ' ' + speciesCoeffs[ispecie] + ')\n')
    f.write(')\n;\n\t\t}\n\t}')

def lindReactionParse(line,f,parser,il):
    f.write('\t\treaction\t\t\"')
    reactionOrder = 0
    reactSpe = line[0].split(parser)
    lhs = reactSpe[0]
    rhs = reactSpe[1]
    lhs = lhs.strip('(+M)')
    rhs = rhs.strip('(+M)')
    reactants = lhs.split('+')
    products = rhs.split('+')
    for j,ireactant in enumerate(reactants):
        if ireactant[0] not in numbers:
            reactionOrder+=1
        else:
            reactionOrder+=int(ireactant[0])
        f.write(ireactant+' ')
        if j != (len(reactants)-1):
            f.write('+ ')
        else:
            pass
    f.write('=')
    for j,iproduct in enumerate(products):
        f.write(' '+iproduct)
        if j != (len(products)-1):
            f.write(' +')
        else:
            pass
    lineLow = lines[il+1]
    lineLow = lineLow.strip()
    lineLow = re.split(r'\s+', lineLow)
    for iLineLow in range(len(lineLow)):
        lineLow[iLineLow] = lineLow[iLineLow].strip('/')
    f.write('";\n')
    f.write('\t\tk0\n\t\t{\n')
    f.write('\t\t\tA\t\t\t' + f'{orderDict[reactionOrder+1]*float(lineLow[1])}' + ';\n')
    f.write('\t\t\tbeta\t\t' + lineLow[2] + ';\n')
    f.write('\t\t\tTa\t\t\t' + f'{float(lineLow[3])*EToTa}' + ';\n\t\t}\n')
    f.write('\t\tkInf\n\t\t{\n')
    f.write('\t\t\tA\t\t\t' + f'{orderDict[reactionOrder]*float(line[1])}' + ';\n')
    f.write('\t\t\tbeta\t\t' + line[2] + ';\n')
    f.write('\t\t\tTa\t\t\t' + f'{float(line[3])*EToTa}' + ';\n\t\t}\n')
    
    f.write('\t\tF\n\t\t{\n' + '\t\t}\n')
    f.write('\t\tthirdBodyEfficiencies\n\t\t{\n')
    f.write('\t\t\tcoeffs\n' + f'{numSpecies}' + '\n(\n')
    line2 = lines[il+2]
    if '/' in line2:
        line2 = line2.strip()
        line2 = re.split(r'\s+', line2)
        for iLine2 in range(len(line2)):
            line2[iLine2] = line2[iLine2].strip('/')
        speciesCoeffs = []
        for ispecie in range(numSpecies):
            speciesCoeffs.append('1')
            for iLine2,modSpec in enumerate(line2):
                if modSpec == species[ispecie]:
                    speciesCoeffs[ispecie] = line2[iLine2+1]
                    break
                else:
                    continue
    else:
        speciesCoeffs = []
        for ispecie in range(numSpecies):
            speciesCoeffs.append('1')
    for ispecie in range(numSpecies):
        f.write('(' + species[ispecie] + ' ' + speciesCoeffs[ispecie] + ')\n')
    f.write(')\n;\n\t\t}\n\t}')

def simpleReactionParse(line,f,parser):
    f.write('\t\treaction\t\t\"')
    reactionOrder = 0
    reactSpe = line[0].split(parser)
    lhs = reactSpe[0]
    rhs = reactSpe[1]
    reactants = lhs.split('+')
    products = rhs.split('+')
    for j,ireactant in enumerate(reactants):
        f.write(ireactant+' ')
        if ireactant[0] not in numbers:
            reactionOrder+=1
        else:
            reactionOrder+=int(ireactant[0])
        if j!= (len(reactants)-1):
            f.write('+ ')
        else:
            continue
    f.write('=')
    for j,iproduct in enumerate(products):
        f.write(' '+iproduct)
        if j!= (len(products)-1):
            f.write(' +')
        else:
            continue
    f.write('";\n')
    f.write('\t\tA\t\t\t\t' + f'{orderDict[reactionOrder]*float(line[1])}' + ';\n')
    f.write('\t\tbeta\t\t\t' + line[2] + ';\n')
    f.write('\t\tTa\t\t\t\t' + f'{float(line[3])*EToTa}' + ';\n\t}')

with open('chem.inp') as fin, open ('reactions','a+') as f:
    for line in fin:
        if 'REACTIONS' in line:
            break
    lines = fin.readlines()
    for il,line in enumerate(lines):
        line = line.strip()
        if 'END' in line:
            break
        else:
            if line[0] == '!' or 'DUPLICATE' in line:
                continue
            elif '+M=>' in line:
                # irreversiblethirdBodyArrheniusReaction
                line = re.split(r'\s+', line)
                f.write('\n\tun-named-reaction-' + f'{numReactions}')
                numReactions += 1
                f.write('\n\t{\n\t\ttype\t\t\tirreversiblethirdBodyArrheniusReaction;\n')
                thirdBodyArrheniusParse(line,f,'=>',il)
            elif '+M=' in line or '+M<=>' in line:
                # reversiblethirdBodyArrheniusReaction
                line = re.split(r'\s+', line)
                f.write('\n\tun-named-reaction-' + f'{numReactions}')
                numReactions += 1
                f.write('\n\t{\n\t\ttype\t\t\treversiblethirdBodyArrheniusReaction;\n')
                if '<=>' in line:
                    thirdBodyArrheniusParse(line,f,'<=>',il)
                else:
                    thirdBodyArrheniusParse(line,f,'=',il)
            elif '(+M)' in line:
                if 'TROE' in lines[il+2]:
                    # reversibleArrheniusTroeFallOffReaction
                    line = re.split(r'\s+', line)
                    f.write('\n\tun-named-reaction-' + f'{numReactions}')
                    numReactions += 1
                    f.write('\n\t{\n\t\ttype\t\t\treversibleArrheniusTroeFallOffReaction;\n')
                    if '<=>' in line:
                        troeReactionParse(line,f,'<=>',il)
                    else:
                        troeReactionParse(line,f,'=',il)
                elif 'SRI' in lines[il+2]:
                    # reversibleArrheniusSRIFallOffReaction
                    line = re.split(r'\s+', line)
                    f.write('\n\tun-named-reaction-' + f'{numReactions}')
                    numReactions += 1
                    f.write('\n\t{\n\t\ttype\t\t\treversibleArrheniusSRIFallOffReaction;\n')
                    if '<=>' in line:
                        sriReactionParse(line,f,'<=>',il)
                    else:
                        sriReactionParse(line,f,'=',il)
                else:
                    # reversibleArrheniusLindemannFallOffReaction
                    line = re.split(r'\s+', line)
                    f.write('\n\tun-named-reaction-' + f'{numReactions}')
                    numReactions += 1
                    f.write('\n\t{\n\t\ttype\t\t\treversibleArrheniusLindemannFallOffReaction;\n')
                    if '<=>' in line:
                        lindReactionParse(line,f,'<=>',il)
                    else:
                        lindReactionParse(line,f,'=',il)
            elif '/' in line:
                continue
            elif '<=>' in line:
                line = re.split(r'\s+', line)
                f.write('\n\tun-named-reaction-' + f'{numReactions}')
                numReactions += 1
                f.write('\n\t{\n\t\ttype\t\t\treversibleArrheniusReaction;\n')
                simpleReactionParse(line,f,'<=>')
            elif '=>' in line:
                line = re.split(r'\s+', line)
                f.write('\n\tun-named-reaction-' + f'{numReactions}')
                numReactions += 1
                f.write('\n\t{\n\t\ttype\t\t\tirreversibleArrheniusReaction;\n')
                simpleReactionParse(line,f,'=>')
            else:
                line = re.split(r'\s+', line)
                f.write('\n\tun-named-reaction-' + f'{numReactions}')
                numReactions += 1
                f.write('\n\t{\n\t\ttype\t\t\treversibleArrheniusReaction;\n')
                simpleReactionParse(line,f,'=')
    f.write('\n}')




# Thermodynamics
weightDict = {'C ':12.011, 'C':12.011, 'H':1.0079, 'H ':1.0079, 'N':14.00674, 'N ':14.00674, 'O':15.999, 'O ':15.999, 'HE':4.0026, 'AR':39.948}
with open('therm.dat') as f, open('thermodynamics','w+') as w:
    while True:
        specie = ''
        coeffs = []
        elements = []
        elemNum = []
        for line in f:
            if line == '\n' or '!' in line:
                continue
            line = line.strip()
            if line[-1] == '1':
                elements.append(line[24:26])
                elemNum.append(int(float(line[26:29])))
                if line[29] != ' ':
                    elements.append(line[29:31])
                    elemNum.append(int(float(line[31:34])))
                if line[34] != ' ':
                    elements.append(line[34:36])
                    elemNum.append(int(float(line[36:39])))
                molWeight = 0.0
                for i,elem in enumerate(elements):
                    molWeight+=elemNum[i]*weightDict[elem]
                line = re.split(r'\s+', line)
                length = len(line)
                specie = line[0]
                Tcommon = line[-2]
                Thigh = line[-3]
                Tlow = line[-4]
                break
        for line in f:
            if line == '\n' or '!' in line:
                continue
            if line[-2] != '4':
                coeffs.append(line[0:15])
                coeffs.append(line[15:30])
                coeffs.append(line[30:45])
                coeffs.append(line[45:60])
                coeffs.append(line[60:75])
            else:   
                coeffs.append(line[0:15])
                coeffs.append(line[15:30])
                coeffs.append(line[30:45])
                coeffs.append(line[45:60])
                break
        if specie == '':
            break
        w.write(specie+'\n{\n\tspecie\n\t{\n\t\tmolWeight\t\t'+f'{molWeight}'+';\n\t}\n')
        w.write('\tthermodynamics\n\t{\n\t\tTlow\t\t\t'+ \
            f'{Tlow};\n\t\tThigh\t\t\t{Thigh};\n\t\tTcommon\t\t\t{Tcommon};\n')
        w.write('\t\thighCpCoeffs\t'+f'( {coeffs[0]} {coeffs[1]} {coeffs[2]} {coeffs[3]} {coeffs[4]}'+ \
            f' {coeffs[5]} {coeffs[6]} );\n')
        w.write('\t\tlowCpCoeffs\t\t'+f'( {coeffs[7]} {coeffs[8]} {coeffs[9]} {coeffs[10]} {coeffs[11]}'+ \
            f' {coeffs[12]} {coeffs[13]} );\n')
        w.write('\t}\n\ttransport\n\t{\n\t\tAs\t\t\t\t1.67212e-06;\n\t\tTs\t\t\t\t170.672;\n\t}\n')
        w.write('\telements\n\t{\n')
        for i,elem in enumerate(elements):
            w.write('\t\t'+elem + '\t\t\t\t' + f'{elemNum[i]};' + '\n')
        w.write('\t}\n}\n')
