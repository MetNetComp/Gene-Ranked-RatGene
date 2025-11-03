import scipy.io as sco
from collections import Counter
import numpy as np
import copy
import os
import argparse
import requests

def splitOutGene(grRule):
    if len(grRule) == 0:
        return None
    else:
        grRule = grRule[0]
    grRule = grRule.replace("(", "")
    grRule = grRule.replace(")", "")
    grRule = grRule.replace("and", "")
    grRule = grRule.replace("or", "")
    glist = grRule.split("  ")
    grlist = Counter(glist)
    return grlist


def addFrequency(frequncylist, grlist):
    for key, value in grlist.items():
        if key in frequncylist:
            frequncylist[key] += grlist[key]
        else:
            frequncylist[key] = grlist[key]
    return frequncylist

def addPropotion(propotionlist, grlist):
    for key in grlist.keys():
        if key in propotionlist:
            propotionlist[key] += 1
        else:
            propotionlist[key] = 1
    return propotionlist

def countTwo(model):
    frequncylist = dict()
    propotionlist = dict()
    nGr = 0
    for i in range(model.grRules.shape[0]):
        grlist = splitOutGene(model.grRules[i][0])
        if not grlist:
            continue
        else:
            nGr += 1
            frequncylist = addFrequency(frequncylist, grlist)
            propotionlist = addPropotion(propotionlist, grlist)
    sfrequncylist = dict(sorted(frequncylist.items(), key=lambda item: item[1], reverse=True))
    spropotionlist = dict(sorted(propotionlist.items(), key=lambda item: item[1], reverse=True))
    return nGr, sfrequncylist, spropotionlist

def maxParseGPR(model):
    # data scale
    nRxn = model.rxns.shape[0]
    nGen = model.genes.shape[0]
    nAux = 1
    parInfo = [0] * nRxn
    nRelation = 0
    nEqual = 0
    gprLabel = np.ones(nRxn)

    # get GPR rules
    for i in range(nRxn):
        gpr = model.grRules[i, 0] # string 
        if  gpr.size <= 0: # empty
            parInfoUnit = np.full((1, 3), "", dtype=object)
            gprLabel[i] = 0
            continue
        gpr = gpr[0]

        if ' ' in gpr: # multiple
            gprStr = "(" + gpr + ")"
            maxUnit = gprStr.count("(")
            nRelation = nRelation + 2 * maxUnit
            parInfoUnit = np.empty((maxUnit, 3), dtype=object)
            unitLabel = 1
            while unitLabel <= maxUnit:
                lbracket = [j for j in range(len(gprStr)) if gprStr[j] == "("]
                rbracket = [j for j in range(len(gprStr)) if gprStr[j] == ")"]
                rpoint = rbracket[0]
                for id in lbracket:
                    if id < rpoint:
                        lpoint = id
                    else:
                        break
                gprUnit = gprStr[lpoint+1:rpoint]
                gprStr = gprStr[0:lpoint] + "aux" + str(nAux) + gprStr[rpoint+1:]

                # assign one to one gene
                if unitLabel == maxUnit:
                    parInfoUnit[unitLabel - 1, 0] = "real" + str(i + 1)
                else:
                    parInfoUnit[unitLabel - 1, 0] = "aux" + str(nAux)
                    nAux += 1
                
                # single and/or boolean
                gprUnitElement = gprUnit.split(" ")
                parInfoUnit[unitLabel - 1, 1] = gprUnitElement[1]
                parInfoUnit[unitLabel - 1, 2] = gprUnitElement[0::2]
                unitLabel += 1
        else: # single
            parInfoUnit = np.empty((1, 3), dtype=object)
            parInfoUnit[0, 0] = "real" + str(i + 1)
            parInfoUnit[0, 1] = ""
            parInfoUnit[0, 2] = gpr
            nEqual += 1
        parInfo[i] = parInfoUnit
    nAux -= 1
    return parInfo, nRxn, nGen, nAux, nRelation, nEqual, gprLabel

def unitScore(itemUnit, score=10):
    geneScore = dict()
    # three types
    if itemUnit[1] == "":
        geneScore[str(itemUnit[2])] = score
    elif itemUnit[1] == "or":
        for gene in itemUnit[2]:
            geneScore[gene] = score
    elif itemUnit[1] == "and":
        for gene in itemUnit[2]:
            geneScore[gene] = score/len(itemUnit[2])
    return geneScore

def findIndex(parInfoUnit, name):
    i = -1
    for i in range(len(parInfoUnit)):
        if parInfoUnit[i][0] == name:
            return i
    return i


def grScore(parInfoUnit, score=10):
    if isinstance(parInfoUnit, int):
        return None
    nRelation = len(parInfoUnit)
    itemUnit = parInfoUnit[nRelation - 1]
    # single gene or single AND/OR relation
    if nRelation == 1:
        geneScore = unitScore(itemUnit, score)
        return geneScore
    # multiple relation
    geneScore = unitScore(itemUnit, score)
    stopLabel = True
    while stopLabel:
        changeLable = 1
        keylist = copy.deepcopy(tuple(geneScore.keys()))
        for key in keylist:
            if key[:3] == "aux":
                changeLable = 0
                id = findIndex(parInfoUnit, key)
                subScore = unitScore(parInfoUnit[id], geneScore[key])
                del geneScore[key]
                for k, v in subScore.items():
                    if k in geneScore:
                        geneScore[k] += v
                    else:
                        geneScore[k] = v
        if changeLable:
            stopLabel = False
    return geneScore

def rxnScore(model, strategy="degree"):
    score = list()
    # naive strategy sum of degrees
    if strategy == "degree":
        for i in range(model.rxns.shape[0]):
            degree = int(10 * sum(model.S[:, i] != 0))
            score.append(degree)
    # count for reversible rxns
    elif strategy == "revdegree":
        for i in range(model.rxns.shape[0]):
            rev = model.rev[i][0] + 1
            multiple_base = rev * 10
            degree = int(multiple_base * sum(model.S[:, i] != 0))
            score.append(degree)
    # count for flux gap
    elif strategy == "flux":
        for i in range(model.rxns.shape[0]):
            flux_gap = model.ub[i][0] - model.lb[i][0]
            degree = np.log(flux_gap + 1) * 10 + 10
            score.append(degree)
    elif strategy == "revflux":
        for i in range(model.rxns.shape[0]):
            rev = model.rev[i][0] + 1
            flux_gap = model.ub[i][0] - model.lb[i][0]
            degree = rev * (np.log(flux_gap + 1) * 10 + 10)
            score.append(degree)
    elif strategy == "fluxdegree":
        for i in range(model.rxns.shape[0]):
            flux_gap = model.ub[i][0] - model.lb[i][0]
            multiple_base = np.log(flux_gap + 1) * 10 + 10
            degree = int(multiple_base * sum(model.S[:, i] != 0))
            score.append(degree)
    elif strategy == "revfluxdegree":
        for i in range(model.rxns.shape[0]):
            rev = model.rev[i][0] + 1
            flux_gap = model.ub[i][0] - model.lb[i][0]
            multiple_base = rev * (np.log(flux_gap + 1) * 10 + 10)
            degree = int(multiple_base * sum(model.S[:, i] != 0))
            score.append(degree)
    else:
        raise ValueError("No such strategy for weighted rxns")
    return score

def countScore(model):
    parInfo, _, _, _, _, _, _ = maxParseGPR(model)
    geneScore = dict()
    for i in range(len(parInfo)):
        gdict = grScore(parInfo[i])
        if not gdict:
            continue
        for key, value in gdict.items():
            if key in geneScore:
                geneScore[key] += value
            else:
                geneScore[key] = value
    geneScore = dict(sorted(geneScore.items(), key=lambda item: item[1], reverse=True))
    return geneScore

def countWeightedScore(model, strategy="degree"):
    parInfo, _, _, _, _, _, _ = maxParseGPR(model)
    geneScore = dict()
    score = rxnScore(model, strategy)
    for i in range(len(parInfo)):
        gdict = grScore(parInfo[i], score[i])
        if not gdict:
            continue
        for key, value in gdict.items():
            if key in geneScore:
                geneScore[key] += value
            else:
                geneScore[key] = value
    geneScore = dict(sorted(geneScore.items(), key=lambda item: item[1], reverse=True))
    return geneScore


def save_countTwo(model, modelname):
    if os.path.exists(f"{modelname}_multiplicity.txt") and os.path.exists(f"{modelname}_frequency.txt"):
        return None, None, None
    nGr, frequncylist, propotionlist = countTwo(model)
    with open(f"{modelname}_multiplicity.txt", "w") as f:
        f.write(f"total grRules:{nGr}\n")
        for key, value in propotionlist.items():
            print(f"{key}:{value/nGr*100}")
            f.write(f"{key}:{value/nGr*100}\n")
    
    with open(f"{modelname}_frequency.txt", "w") as f:
        for key, value in propotionlist.items():
            f.write(f"{key}:{value}\n")
            print(f"{key}:{value}")
    return nGr, frequncylist, propotionlist

def save_score(model, modelname):
    if os.path.exists(f"{modelname}_logic.txt"):
        return None
    geneScore = countScore(model)
    with open(f"{modelname}_logic.txt", "w") as f:
        for key, value in geneScore.items():
            f.write(f"{key}:{value}\n")
            print(f"{key}:{value}")
    return geneScore

def save_weightedScore(model, modelname, strategy="degree"):
    if os.path.exists(f"{modelname}_score_{strategy}.txt"):
        return None
    geneScore = countWeightedScore(model, strategy)
    with open(f"{modelname}_score_{strategy}.txt", "w") as f:
        for key, value in geneScore.items():
            f.write(f"{key}:{value}\n")
            print(f"{key}:{value}")
    return geneScore

def main(modelname, strategy):
    readm = sco.loadmat(f'data/{modelname}.mat', struct_as_record=False)
    model = readm[modelname][0][0]
    if strategy == "multiplicity" or strategy == "frequency" or strategy == "all":
        _, _, _ = save_countTwo(model, modelname)
        print(f"saved in {modelname}_{strategy}.txt")

    if strategy == "logic" or strategy == "all":
        _ = save_score(model, modelname)
        print(f"saved in {modelname}_logic.txt")

    scorelsit = ["degree", "revdegree", "flux", "revflux", "fluxdegree", "revfluxdegree"]   
    if strategy in scorelsit:
        _ = save_weightedScore(model, modelname, strategy)
        print(f"saved in {modelname}_score_{strategy}.txt")
    
    if strategy == "all":
        for s in scorelsit:
            _ = save_weightedScore(model, modelname, s)
            print(f"saved in {modelname}_score_{s}.txt")
    print("end")


if __name__ == "__main__":     
    # parse params
    parser = argparse.ArgumentParser(description='Hyper Parameters for Gene-Ranked RatGene')
    parser.add_argument(
        '--model', 
        type=str, 
        help='name of the .mat model file downloaded from BiGG, should be inside of ./data/',
        default='e_coli_core'
        )
    parser.add_argument(
        '--strategy', 
        type=str, 
        help='name of the strategy used for ranking the gene imporatance',
        default='frequency'
        )

    args = parser.parse_args()
    modelName = args.model
    strategy = args.strategy

    strategies = ["multiplicity", "frequency", "logic", "degree", "revdegree", "flux", "revflux", "fluxdegree", "revfluxdegree", "all"]
    if strategy not in strategies:
        raise ValueError(f"{strategy} name is not correct.")
    
    # load .mat data
    modelName = "iJR904"
    dataFolder = os.path.join(os.getcwd(), "data")
    modelNames = os.listdir(os.path.join(dataFolder))
    modelChoice = [f.split(".")[0] for f in modelNames]
    if modelName not in modelChoice:
        url = f"http://bigg.ucsd.edu/static/models/{modelName}.mat"
        response = requests.get(url)
        # download dataset
        if response.status_code == 200:
            with open(f"data/{modelName}.mat", "wb") as f:
                f.write(response.content)
            print(f"successfully downloaded {modelName}.mat")
        else:
            raise ValueError(f"downlaod error, status: {response.status_code}, please download the dataset from BiGG manually.")
    
    # Gene-Ranked RatGene
    main(modelName, strategy)
