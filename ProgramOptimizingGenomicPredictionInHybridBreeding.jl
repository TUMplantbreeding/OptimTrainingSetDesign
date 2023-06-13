using Dates
using Random
using Distributions
using Printf
using DataFrames
using LinearAlgebra
using Plots
using CSV
Random.seed!(31415)

# Definitions of the parameters
nRepPar=20
nRepTrait=50
nTested = [  #c and nTS in the first respectively in the second column.
    1 108; 1 96; 1 84; 1 72; 1 60; 1 48; 1 36; 1 24; 1 12;
    2 108; 2 96; 2 84; 2 72; 2 60; 2 48; 2 36; 2 24; 2 12;
    4 108; 4 96; 4 84; 4 72; 4 60; 4 48; 4 36; 4 24; 4 12;
];

genoSzen =[ # nQTL, nD0, sdDom, meanDom, targetRatioSCAG;
    1000 500 .33  1  .22;
    1000   0 .15 .5  .06;
     300   0 .06 .18 .01
]; # in the columns are
# column 1 number of QTL
# column 2 number of QTL with only dominance effects
# column 3 standard deviation of dominance effects
# column 4 mean of standard deviation of dominance effects
# column 5 target ratio of SCA to G 

#function definitions

function rhoPi(id, L, G) 
# approximation of expected accuracy for individuals in set id, which is a vector with indices of the effects in the set, 
# according to Eqs. (13) and (B6)
# L = BVB' = Var(uHat) given in Eq. (6)
# G = Var(u) given in Eq. (3), where u is the vector of random effects in the model
    n = size(id,1)
    trNum = 0
    for i in id
        trNum += sum(L[i,i])
    end
    trNum -= sum(L[id,id])/n
    trDenom = 0
    for i in id
        trDenom += G[i, i]
    end
    trDenom -= sum(G[id,id])/n
    sqrt(trNum/trDenom)
end #rhoPi


function rho_Pi(L, G, idM0, idM1, idF0, idF1, idH0, idH1_M, idH1_F, idH2, idH3, idH)
# Calling function rhoPi(id, L, G) for 10 different sets 
    vec([rhoPi(idM0, L, G);
            rhoPi(idM1, L, G);
            rhoPi(idF0, L, G);
            rhoPi(idF1, L, G);
            rhoPi(idH0, L, G);
            rhoPi(idH1_M, L, G);
            rhoPi(idH1_F, L, G);
            rhoPi(idH2, L, G);
            rhoPi(idH3, L, G);
            rhoPi(idH, L, G)])
end # rho_Pi for 10 sets

function rho_Pi(L, G, idH0, idH1_M, idH1_F, idH2, idH3, idH)
# Calling function rhoPi(id, L, G) for 6 different sets     
    vec([rhoPi(idH0, L, G);
         rhoPi(idH1_M, L, G);
         rhoPi(idH1_F, L, G);
         rhoPi(idH2, L, G);
         rhoPi(idH3, L, G);
         rhoPi(idH, L, G)])
end # rho_Pi for 6 sets

function rho_Pi(L, G, idM0, idM1, idF0, idF1)
# Calling function rhoPi(id, L, G) for 4 different sets     
    vec([rhoPi(idM0, L, G);
            rhoPi(idM1, L, G);
            rhoPi(idF0, L, G);
            rhoPi(idF1, L, G)])
end # rho_Pi for 4 sets

function rhoPiW(id, L, G)
# approximation of expected accuracy for hybrids in set id, which is a vector with indices of the hybrids in the set, 
# according to Eq. (14) 
# L = BVB' = Var(uHat) given in Eq. (6)
# G = Var(u) given in Eq. (3), where u is the vector of random effects in the model
    n = size(id,1)
    trNum = 0
    for i in id
        trNum += sum(L[[i ped[i,1] ped[i,2]], [i ped[i,1] ped[i,2]]])
    end
    trNum -= sum(L[[id ped[id,1] ped[id,2]], [id ped[id,1] ped[id,2]]])/n
    trDenom = 0
    for i in id
        trDenom += G[i, i]
    end
    for i in ped[id,1]
        trDenom += G[i, i]
    end
    for i in ped[id,2]
        trDenom += G[i, i]
    end
    trDenom -= sum(G[[id ped[id,1] ped[id,2]], [id ped[id,1] ped[id,2]]])/n
    sqrt(trNum/trDenom)
end # rhoPiW

function rho_PiW(L, G, idH0, idH1_M, idH1_F, idH2, idH3, idH)
# Calling function rhoPiW(id, L, G) for 6 different sets     
    vec([rhoPiW(idH0, L, G);
         rhoPiW(idH1_M, L, G);
         rhoPiW(idH1_F, L, G);
         rhoPiW(idH2, L, G);
         rhoPiW(idH3, L, G);
         rhoPiW(idH, L, G)])
end # rho_PiW for different sets of hybrid values

function pedFill(idM, idF)
# construct a matrix with the pedigree.
# parent1 in the first column and parent2 in the second column
# in the first rows are zeros for the inbred lines, so that hybrid identifier can used directly
    ped = Array{Int, 2}(undef,(size(idM,1)+1)*(size(idF,1)+1)-1, 2)
    row = size(idM,1)+size(idF,1)
    for i in idM, j in idF
        row += 1
        ped[row,1] = i
        ped[row,2] = j
    end
    return ped
end #pedFill

function corForIds(vec1, vec2, idM0, idM1, idF0, idF1, idH0, idH1_M, idH1_F, idH2, idH3, idH)
#? What eq or table in paper?    
# calculate correlations between members of vec1 and vec2
# for 10 different sets
    corVals= zeros(10)
    corVals[1] = cor(vec1[idM0], vec2[idM0])
    corVals[2] = cor(vec1[idM1], vec2[idM1])
    corVals[3] = cor(vec1[idF0], vec2[idF0])
    corVals[4] = cor(vec1[idF1], vec2[idF1])
    corVals[5] = cor(vec1[idH0], vec2[idH0])
    corVals[6] = cor(vec1[idH1_M], vec2[idH1_M])
    corVals[7] = cor(vec1[idH1_F], vec2[idH1_F])
    corVals[8] = cor(vec1[idH2], vec2[idH2])
    corVals[9] = cor(vec1[idH3], vec2[idH3])
    corVals[10] = cor(vec1[idH], vec2[idH])
    return vec(corVals)
end # corForIds for 10 id vectors

function corForIds(vec1, vec2, idH0, idH1_M, idH1_F, idH2, idH3, idH)
# calculate correlations between members of vec1 and vec2
# for 6 different sets
    corVals= zeros(6)
    corVals[1] = cor(vec1[idH0], vec2[idH0])
    corVals[2] = cor(vec1[idH1_M], vec2[idH1_M])
    corVals[3] = cor(vec1[idH1_F], vec2[idH1_F])
    corVals[4] = cor(vec1[idH2], vec2[idH2])
    corVals[5] = cor(vec1[idH3], vec2[idH3])
    corVals[6] = cor(vec1[idH], vec2[idH])
    return vec(corVals)
end # corForIds for 6 id vectors

function corForIds(vec1, vec2, idM0, idM1, idF0, idF1)
# calculate correlations between members of vec1 and vec2
# for 4 different sets
    corVals= zeros(4)
    corVals[1] = cor(vec1[idM0], vec2[idM0])
    corVals[2] = cor(vec1[idM1], vec2[idM1])
    corVals[3] = cor(vec1[idF0], vec2[idF0])
    corVals[4] = cor(vec1[idF1], vec2[idF1])
    return vec(corVals)
end  # corForIds for 4 id vectors

function getIndices(idM0, idM1, idF0, idF1, nSel, nC)
# get indices for the different hybrid groups
# indices are needed for calculation of corelations and accuracies
    nH3 = nSel*nC # get indices for H3
    pedH3 = zeros(Int64, nH3, 2) # pedigree for hybrids in training data (H3)
    idH3  = zeros(Int64, nH3)    # index vector for hybrids in training data (H3)
    k = 0
    for i=1:nSel
        ii = i
        for j=1:nC
            jj = ((i-1) + (j-1))%nSel + 1
            k += 1
            pedH3[k,1] = idM1[ii]
            pedH3[k,2] = idF1[jj]
            idH3[k] = nM + nF  + (idM1[ii]-1)*nF + (idF1[jj] - nM) # index of hybrid from male idMsel[ii] 
                                                                   # and female idFsel[ii]
                                                                   # hybrid indices start from nM+nF+1
                                                                   # female indices start from nM+1
        end
    end
    # get indices for H0
    nH0 = length(idM0)*length(idF0)
    idH0 = zeros(Int64, nH0)
    k = 0
    for i in idM0
        for j in idF0
          k+=1
          idH0[k]  = nM + nF  + (i-1)*nF    +   (j - nM)
        end
    end
    # get indices for H1_F
    nH1_F = length(idM0)*length(idF1)
    idH1_F = zeros(Int64, nH1_F)
    k = 0
    for i in idM0
        for j in idF1
          k+=1
          idH1_F[k] = nM + nF  + (i-1)*nF + (j - nM)
        end
    end
    # get indices for H1_M
    nH1_M = length(idM1)*length(idF0)
    idH1_M = zeros(Int64, nH1_M)
    k = 0
    for i in idM1
        for j in idF0
          k+=1
          idH1_M[k] = nM + nF  + (i-1)*nF + (j - nM)
        end
    end
    # get indices for H2
    nH2 = length(idM1)*length(idF1)
    idH2 = zeros(Int64, nH2)
    k = 0
    for i in idM1
        for j in idF1
          k+=1
          idH2[k] = nM + nF + (i-1)*nF + (j - nM)
        end
    end
    idH2 = setdiff(idH2,idH3);
   return(idH0, idH1_M, idH1_F, idH2, idH3, pedH3)
end # getIndices
#end of function definitions

# script for reading in all data

Random.seed!(31415)
# definitions for parental lines
nM  = 145 # number of males !!! In this case these are Dent lines, data set 1!!!
nF  = 111 # number of females !!! In this case these are Flint lines, data set 1 !!!

idM = collect(1:nM)                 # index vector for males
idF = collect(nM + 1:(nM+nF))       # index vector for females
idH = collect(nM+nF+1:nM+nF+nM*nF); # index vector for hybrids
# pedigree for all possible hybrids 
ped = pedFill(idM, idF)

df = CSV.read("InputGeno.csv", DataFrame; header=0)
genMatInput = Matrix{Int64}(df)
dfPos = CSV.read("map_pruned.csv", DataFrame; header=1)

genMatFull = zeros(nM+nF+nM*nF,size(genMatInput,2))
genMatFull[1:nM+nF, :] = genMatInput*2

# generate genotypes of all possible hybrids
for i in idH
    male   = ped[i,1]
    female = ped[i,2]
    genMatFull[i,:] = genMatInput[male,:] + genMatInput[female,:]
end

#split genotypic matrices in Marker for QTL and SNPs for kinship this step is not done for data set 2
QTLpos = zeros(Int64, 3000)
chr = zeros(Int64, 10)
for i in 1:10
    chr[i] = round(Int64, 3000/size(genMatInput,2)*sum(dfPos[:,2] .== i-1))
    QTLpos[sum(chr[1:i])+1:sum(chr[1:i])+round(Int64, 3000/size(genMatInput,2)*sum(dfPos[:,2] .== i))] = 
            sample((1:size(genMatInput,2))[dfPos[:,2] .== i], round(Int64, sum(dfPos[:,2] .== i)*3000/size(genMatInput,2)), replace=false)
end
markMat = genMatFull[:, QTLpos]
nLoci = size(markMat,2)
genMat = genMatFull[:, setdiff(1:size(genMatFull,2),QTLpos)]
# Calculation of kinship matrices. This is like simple matching!

# centering genotypes for males
genM = genMat[idM,:]
meanM = mean(genM,dims=1)
genM .-= repeat(meanM,nM)# mean(genM,dims=1)
#genM = genM[1:end-1,:] # skip last line for full rank

# centering genotypes for females
genF = genMat[idF,:]
meanF = mean(genF,dims=1)
genF .-= repeat(meanF, nF) # mean(genF,dims=1)
#genF = genF[1:end-1,:] # skip last line for full rank

# centering genotypes for hybrids
genH = genMat[idH,:]
meanH = mean(genH,dims=1)
genH .-= repeat(meanH, nM*nF) # mean(genH,dims=1)

# computing kinship matrices
K_M = genM*genM' ./ (meanM * (2 .-meanM)')  # GCA males
K_F = genF*genF' ./ (meanF * (2 .-meanF)')  # GCA females
K_H = kron(K_M,K_F)               # SCA

K_Mr = K_M[1:end-1, 1:end-1] # get rid of the last row and column because of linear dependencies
K_Fr = K_F[1:end-1, 1:end-1]
K_Hr = kron(K_Mr, K_Fr)      # relationship for the hybrids

iK_M = inv(K_Mr)
iK_F = inv(K_Fr)
iK_H = kron(iK_M, iK_F)
indexKinH = vec(reshape(idH, nF,nM)[1:end-1, 1:end-1])
println("Module input done.\n")

# end of script for reading in all data
##########################################################################################################################
# loop for different trait architecture set to 1 to shorten the program runtime
# if used as loop, be aware to expand the containers for the results and delete the comment before the end statement!
cntSzen = 1
#for cntSzen in 1:size(genoSzen,1) 
    ######################################################################################################################
    # module genetic values in chart Module A in the paper
    println("Counter for Trait samples nRepTrait: ",nRepTrait)
    println(string("genotypic Scenorio: nQTL:", Int(genoSzen[cntSzen, 1]), " nD0: ", 
            Int(genoSzen[cntSzen, 2]), " sdD: ", genoSzen[cntSzen, 3], " muD: ", genoSzen[cntSzen, 4]))
    println(string("                    Target SCA/G: ", genoSzen[cntSzen, 5]))
    nQTL    = Int64(genoSzen[cntSzen, 1])   # number of QTL
    nD0     = Int64(genoSzen[cntSzen, 2])  # number of QTL with only dominance effects
    sdDom   = genoSzen[cntSzen, 3]  # standard deviation for dominance effects
    meanDom = genoSzen[cntSzen, 4]  # mean for dominance effects
    TargetRatioSCAG = genoSzen[cntSzen, 5]
    Random.seed!(31415)
    distrib = Gamma(.4, 5/3)
    gGCASCAvar = zeros(nRepTrait, 4)           # container for GCA and SCA variances
    TGV = zeros(size(markMat,1), nRepTrait)    # container for the true breeding values
    GCASCA = zeros(size(markMat,1), nRepTrait) # container for GCA and SCA values
    PV =  zeros(size(markMat,1), nRepTrait, 2) # container for phenotypic values, one for each hybrid
    meanSCArat = zeros(2) # needed as counter

    for cnt in 1:size(TGV,2)
        ratSCAg = 0; cnt2 = 0
        while ratSCAg < .9*TargetRatioSCAG || ratSCAg > 1.1*TargetRatioSCAG
            if cnt2 > 50 # this is needed, if you change parameters and your target ratio does not fit your expectation.
                println(string("Mean of ratio SCA/G: ", meanSCArat[2], " runs: ", meanSCArat[1], "  ", meanSCArat[1]/nRepTrait))
                error()
            end
            cnt2 += 1
            meanSCArat[1] +=1
            idQTL  = sample(1:nLoci, nQTL, replace=false) # position of QTL
            idD0   = sample(1:nQTL, nD0, replace=false)   # position of QTL with only dominance
            gammaEff = rand(distrib, nQTL)                # SNP effects
            idAdd  = idQTL[Not(idD0)]                     # positions of QTL with additive and dominance effects
            signAdd = 2*rand(Bernoulli(0.5), nQTL) .-1    # a sign is needed, because SNP coding is done by allele frequency!
            addEff = signAdd[Not(idD0)].*gammaEff[Not(idD0)] # actual additive effects

            domEff = gammaEff
            k = randn(nQTL) * sdDom .+ meanDom            # degree of dominance
            for i in 1:nQTL
            domEff[i,1] *= k[i,1]                         # actual dominance effects
            end

            TGV[:,cnt] = markMat[:,idAdd] * addEff .+ (markMat[:,idQTL] .== 1.0) * domEff # true genetic values for all members of the data set

        # calculate GCA and SCA and overall
            genH = reshape(TGV[idH,cnt],nF,nM)
            meanH = mean(genH)
            genH .-= meanH
            GCASCA[idM,cnt] = gcaM = vec(mean(genH,dims=1)) # GCA values for the male lines
            GCASCA[idF,cnt] = gcaF = vec(mean(genH,dims=2)) # GCA values for the female lines
            sca = genH .- ones(nF,1)*gcaM'                  # subtract from the hybrid values gcaM
            sca .-=  gcaF * ones(1,nM)                      # and subtract from that gcaF
            GCASCA[idH,cnt] = vec(sca);                     # SCA values for the hybrids
            # GCA and SCA variances, included covariances, reduced to avoid dependencies
            gGCASCAvar[cnt, 1] = GCASCA[idM[1:end-1],cnt]'iK_M*GCASCA[idM[1:end-1],cnt]/(size(gcaM,1)-1) # variance of GCA male
            gGCASCAvar[cnt, 2] = GCASCA[idF[1:end-1],cnt]'iK_F*GCASCA[idF[1:end-1],cnt]/(size(gcaF,1)-1) # variance of GCA female
            gGCASCAvar[cnt, 3] = GCASCA[indexKinH,cnt]'iK_H*GCASCA[indexKinH,cnt]/size(indexKinH,1)      # variance of SCA
            gGCASCAvar[cnt, 4] = var(TGV[idH,cnt])                                                       # variance of TBV without covariances
            ratSCAg = gGCASCAvar[cnt, 3]/(gGCASCAvar[cnt, 1]+gGCASCAvar[cnt, 2]+gGCASCAvar[cnt, 3])      # ratio of SCA to total genetic variance
            meanSCArat[2] += (ratSCAg - meanSCArat[2]) / meanSCArat[1]                                   # calulate the mean of all runs
        end
        PV[:,cnt, 1] = TGV[:,cnt] + randn(size(TGV[:,cnt],1))*sqrt(gGCASCAvar[cnt, 4])*.25; # phenotypes with heritability 0.8
        PV[:,cnt, 2] = TGV[:,cnt] + randn(size(TGV[:,cnt],1))*sqrt(gGCASCAvar[cnt, 4])*1.5; # phenotypes with heritability 0.4
    end
    println(string("Mean of ratio SCA/G: ", meanSCArat[2], " runs: ", meanSCArat[1], "  ", meanSCArat[1]/nRepTrait))
    println("Module genetic values done.\n")

    # script for parent sampling, module B in the chart
    println("Counter for Parents samples nRepPar: ",nRepPar)
    Random.seed!(31415)
    idMsamp = Array{Int, 2}(undef,nM, nRepPar)
    idFsamp = Array{Int, 2}(undef,nF, nRepPar)
    for cnt in 1:nRepPar  # get parental lines for the hybrids going to be tested
        idMsamp[:,cnt] = sample(idM,nM,replace=false) # sample of nTS males
        idFsamp[:,cnt] = sample(idF,nF,replace=false) # sample of nTS females
    end
    println("Module parents done.\n")

    # module C, calulations of predictions
    println("Module evaluation")
    popAcc = zeros(size(TGV,2), size(idMsamp,2),size(nTested,1),2,28) # container for r_hat
    predAcc = zeros(size(TGV,2), size(idMsamp,2),size(nTested,1),2,28) # container for r_a
    eqn16 = zeros(size(TGV,2), size(idMsamp,2),size(nTested,1),2,10) # container for r_tilda

    for cntnTSnC in 1:size(nTested,1) # loop for different nTS and number of crosses
        nTS = nTested[cntnTSnC, 2] # number of selected parents
        nC  = nTested[cntnTSnC, 1] # number of crosses per parent
        for cntH2 in 1:2 # loop for different heritabilities
            nH3 = nTS*nC # number of hybrids in training 
            for cntPar in 1:nRepPar, cntTrait in 1:nRepTrait #loop for samples of parents and traits
                println("Loop: nTS", nTS, " nC ", nC, " Par: ", cntPar, " trait: ", cntTrait, " H2: ", [.8 .4][cntH2], " time: ", Dates.now())
                flush(stdout)

                # get parental lines for the hybrids going to be tested
                idM1 = idMsamp[1:nTS,cntPar] # sample of nTS males, for the specific loop
                idF1 = idFsamp[1:nTS,cntPar] # sample of nTS females, for the specific loop
                idM0 = setdiff(idM,idM1)
                idF0 = setdiff(idF,idF1)
                idH0, idH1_M, idH1_F, idH2, idH3, pedH3 = getIndices(idM0, idM1, idF0, idF1, nTS, nC)
                # spilt the hybrids in the different groups

                # construct G for this model (equation 3)
                G = zeros(nM+nF+nM*nF,nM+nF+nM*nF)
                G[idM,idM] = K_M*gGCASCAvar[cntTrait, 1] # sigma2gcaM 
                G[idF,idF] = K_F*gGCASCAvar[cntTrait, 2] # sigma2gcaF
                G[idH,idH] = K_H*gGCASCAvar[cntTrait, 3] # sigma2sca;
            

                # get Z for model with gcaM, gcaF and sca (equation 3)
                Z = zeros(nH3,nM+nF+nM*nF)
                for i=1:nH3
                    m  = pedH3[i,1]
                    f  = pedH3[i,2] 
                    ii = idH3[i]  
                    Z[i,m]  = 1.0
                    Z[i,f]  = 1.0
                    Z[i,ii] = 1.0
                end
                V = Z*G*Z' + I*gGCASCAvar[cntTrait,4] * [.25 1.5][cntH2] # equation 2-3
                y = PV[idH3, cntTrait, cntH2]; # phenotype vector
                X = ones(nH3,1) # mean vector
                Vi = inv(V)
                # defined for equation 4
                B = G*Z'(Vi*(I - X*XVX*X'Vi))
                # gHat in equation 5; in the paper ghat and shat for GCA and SCA are given seperately
                gHat = B*y;

                # get W for all hybrids to calculate the hybrid values (equation 8)
                W = zeros(nM*nF,nM+nF+nM*nF)
                for i in idH
                    ii = i - nM - nF
                    male   = ped[i,1]
                    female = ped[i,2] 
                    W[ii, male  ]  = 1.0
                    W[ii, female]  = 1.0
                    W[ii ,i     ]  = 1.0
                end
                # equation 9
                hHat = W*gHat;
                # L is defined in Appendix B, needed for equations 16 and 17
                L= B*V*B';
                
                # prediction accuracy (equation 11 and 12)
                # columns 1-10 (equation 11), GCA and SCA only, columns 11-16 hybrid values (equation 12)
                # columns 17-28 (equation 11) for GCA values of hybrids
                predAcc[cntTrait, cntPar, cntH2, cntnTSnC, 1:4] = 
                    corForIds(GCASCA[:, cntTrait], gHat, idM0, idM1, idF0, idF1)
                predAcc[cntTrait, cntPar, cntH2, cntnTSnC, 5:10] = 
                    corForIds(GCASCA[:, cntTrait], gHat, idH0, idH1_M, idH1_F, idH2, idH3, idH)
                predAcc[cntTrait, cntPar, cntH2, cntnTSnC, 11:16] = 
                    corForIds(TGV[idH, cntTrait], hHat, idH0.-nM.-nF, idH1_M.-nM.-nF, idH1_F.-nM.-nF, 
                                    idH2.-nM.-nF, idH3.-nM.-nF, idH.-nM.-nF)
                predAcc[cntTrait, cntPar, cntH2, cntnTSnC, 17:22] = 
                    corForIds(GCASCA[:, cntTrait], gHat, ped[idH0,1], ped[idH1_M,1], ped[idH1_F,1], ped[idH2,1], ped[idH3,1], ped[idH,1])
                predAcc[cntTrait, cntPar, cntH2, cntnTSnC, 23:28] = 
                    corForIds(GCASCA[:, cntTrait], gHat, ped[idH0,2], ped[idH1_M,2], ped[idH1_F,2], ped[idH2,2], ped[idH3,2], ped[idH,2])
 
                # r_hat (equations 13 and 14)
                # columns 1-10 (equation 13) GCA and SCA only, columns 11-16 hybrid values (equation 14)
                # columns 17-28 (equation 13) for GCA values of hybrids
                popAcc[cntTrait, cntPar, cntH2, cntnTSnC, 1:10] = rho_Pi(L, G, idM0, idM1, idF0, idF1, idH0, idH1_M, idH1_F, idH2, idH3, idH)
                popAcc[cntTrait, cntPar, cntH2, cntnTSnC, 11:16] = rho_PiW(L, G, idH0, idH1_M, idH1_F, idH2, idH3, idH)
                popAcc[cntTrait, cntPar, cntH2, cntnTSnC, 17:22] = rho_Pi(L, G, ped[idH0,1], ped[idH1_M,1], ped[idH1_F,1], ped[idH2,1], ped[idH3,1], ped[idH,1])
                popAcc[cntTrait, cntPar, cntH2, cntnTSnC, 23:28] = rho_Pi(L, G, ped[idH0,2], ped[idH1_M,2], ped[idH1_F,2], ped[idH2,2], ped[idH3,2], ped[idH,2])

                # equation 16
                tau = vec(gGCASCAvar[cntTrait, 1:3])/sum(gGCASCAvar[cntTrait, 1:3])
                eqn16[cntTrait, cntPar, cntH2, cntnTSnC, 5:10] = sqrt(vec(predAcc[cntTrait, cntPar, cntH2, cntnTSnC, 17:22])^2/tau[1] + 
                                                            vec(predAcc[cntTrait, cntPar, cntH2, cntnTSnC, 23:28])^2/tau[2] + 
                                                            vec(predAcc[cntTrait, cntPar, cntH2, cntnTSnC, 5:10])^2/tau[3])
            end # end of loop for samples of parents and traits
        end # end of loop for different heritabilities
    end # end of loop for different  nTS and number of crosses
#end # the end of the for loop for genetic architecture
println("Done!")
