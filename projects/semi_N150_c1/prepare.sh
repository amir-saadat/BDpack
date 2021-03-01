dataDir="data"
semidiluteSubDir=(Ree Rg CoM dump m rheol st rst)

BDpack=$BDPACKROOT/bin/BDpack

if [[ ! -d "$dataDir" ]]
then
  echo "data directory doesn't exist. Creating now.."
  mkdir ./$dataDir
  echo "data direcotry created"
else
  echo "data directory exists"
fi

if grep -q semidilute_bs "input.dat"
then
  for ((isubdir=0; isubdir<${#semidiluteSubDir[*]}; isubdir++))
  do
    if [[ ! -d $dataDir/${semidiluteSubDir[isubdir]} ]]
    then
      echo "Following subdir doesn't exist. Creating now..:" ${semidiluteSubDir[isubdir]}
      mkdir ./$dataDir/${semidiluteSubDir[isubdir]}
    else
      echo "Following subdir exists:" ${semidiluteSubDir[isubdir]}
    fi
  done
fi
