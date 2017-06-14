# script for getting simple polymorphism and divergence statistics from the element and flank EM input files

set -e
set -o pipefail

elementEMInputFile=$1
flankEMInputFile=$2
outFile=$3
thresList="$4"
setName="$5"

echo $setName

numSamples=`cat $flankEMInputFile | awk '{if($1=="samples"){print $2;}}'`
#echo $numSamples;

# get neutral SFS from the flank file

neutralSFS=`\
   cat $flankEMInputFile | awk -v OFS="\t" -v numSamples=$numSamples '\
      BEGIN{ \
         for(c=1; c<numSamples; c++) { counts[c]=0.0; } \
         lambda = 0.0; theta = 0.0; \
         polyCount = 0.0; \
      } \
      { \
         if($1=="block") { \
            lambda = $6; theta = $4; \
         } else if($1=="site" && $3 == "P") { \
            polyCount++;
            qAmaj = (1-lambda)*$4 + lambda/3.0*(1-$4); \
            qAmin = (1-lambda)*$5 + lambda/3.0*(1-$5); \
            qTot  = qAmaj + qAmin; \
            pAmaj = qAmaj / qTot; \
            pAmin = qAmin / qTot; \
            counts[$7] += pAmaj; \
            counts[$6] += pAmin; \
         } \
      } \
      END{ \
         countLine=""; \
         for(c=1; c<numSamples; c++) { countLine=countLine"\t"(counts[c]/polyCount); } \
         print countLine; \
      }' \
`
#echo $neutralSFS


# parse element file

elementSummary=`\
   cat $elementEMInputFile | awk -v OFS="\t" -v numSamples=$numSamples '\
      BEGIN{ \
         for(c=1; c<numSamples; c++) { counts[c]=0.0; } \
         lambda = 0.0; theta = 0.0; \
         polyCount  = 0.0; \
         nodivCount = 0.0; \
         divCount   = 0.0; \
         polyExp    = 0.0; \
         nodivExp   = 0.0; \
         divExp     = 0.0; \
         wattersonA=0; for(i=1; i<numSamples-1; i++) {wattersonsA+=(1.0/i);} \
      } \
      { \
         if($1=="block") { \
            lambda = $6; theta = $4; \
            ePoly=theta*wattersonsA; \
            eNodiv=(1-ePoly)*(1-lambda); \
            eDiv=(1-ePoly)*lambda; \
         } else if($1=="site") { \
            if($3 == "P" || $3 == "M") { \
               polyExp  += ePoly; \
               nodivExp += eNodiv; \
               divExp   += eDiv; \
               if($3 == "P") { \
                  polyCount++;
                  qAmaj = (1-lambda)*$4 + lambda/3.0*(1-$4); \
                  qAmin = (1-lambda)*$5 + lambda/3.0*(1-$5); \
                  qTot  = qAmaj + qAmin; \
                  if(qTot==0) { print "P qTot=0 , lambda="lambda" , theta="theta" , "$0; } \
                  pAmaj = qAmaj / qTot; \
                  pAmin = qAmin / qTot; \
                  counts[$7] += pAmaj; \
                  counts[$6] += pAmin; \
               } else { \
                  qNodiv = (1-lambda)*$4; \
                  qDiv   = lambda/3.0*(1-$4); \
                  qTot  = qNodiv + qDiv; \
                  if(qTot==0) { print "M qTot=0 , lambda="lambda" , theta="theta" , "$0; } \
                  pDiv = qDiv / qTot; \
                  pNodiv = qNodiv / qTot; \
                  divCount += pDiv; \
                  nodivCount += pNodiv; \
               } \
            } \
         } \
      } \
      END{ \
         countLine=nodivCount"\t"divCount"\t"polyCount; \
         for(c=1; c<numSamples; c++) { countLine=countLine"\t"(counts[c]/polyCount); } \
         countLine=countLine"\tENDLINE\t"nodivExp"\t"divExp"\t"polyExp; \
         print countLine; \
      }' \
`
header="nodiv\tdiv\tpoly"
for i in `seq 1 $[$numSamples-1]`
do
   header=${header}"\tf["$i"]"
done
tableRow1=`echo ${elementSummary}${neutralSFS} | tr "\n" " " | sed -e 's/ ENDLINE /\\n/' | head -n1`
tableRow2=`echo ${elementSummary}${neutralSFS} | tr "\n" " " | sed -e 's/ ENDLINE /\\n/' | tail -n1`

if [ ""$outFile == "" ]
then
   echo -e $header
   echo $tableRow1
   echo $tableRow2
else
   echo -e $header  > $outFile.txt
   echo $tableRow1 >> $outFile.txt   
   echo $tableRow2 >> $outFile.txt
   
   for thres in $thresList
   do
      echo "thres "$thres
      R --vanilla < ./scripts/executeEM/plotPolyDivStats.r --args $outFile.txt $outFile.f${thres}.pdf $thres "$setName" \
      &> /dev/null
   done
fi


