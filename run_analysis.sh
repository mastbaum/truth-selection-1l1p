if [$# == 3]
then
  Dataset=$1
  directory=$2
  identifier=$3
elif [$# == 2]
then
  directory=$1
  identifier=$2
  Dataset=1
elif [$# == 1]
then
  identifier=$1
  directory="."
  Dataset=1
fi

python scripts/selection.py  $directory/nue_${identifier}.root $Dataset "/pnfs/uboone/scratch/users/mastbaum/comparison/nue_MC_out/v06_47_00/19923250_*/prod*.root"
python scripts/selection.py  $directory/inclusive_${identifier}.root $Dataset "/pnfs/uboone/scratch/users/mastbaum/comparison/inclusive_MC_out/v06_47_00/19923235_*/prod*.root"
python scripts/selection.py  $directory/signal_${identifier}.root $Dataset "/pnfs/uboone/scratch/users/mastbaum/comparison/nue_signal_MC_out/v06_47_00/23371063_*/prod*.root"

root -l -x -q "scripts/utils/add.C(\"nue_${identifier}.root\",\"inclusive_${identifier}.root\",\"background_${identifier}.root\")" 

mkdir cov_background_${identifier}
mkdir cov_signal_${identifier}

python scripts/cov.py background_${identifier}.root -d cov_background_${identifier}
python scripts/cov.py signal_${identifier}.root -s -d cov_signal_${identifier}

python scripts/analysis/chi2.py cov_background_${identifier} cov_signal_${identifier}

