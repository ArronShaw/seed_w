# @Author: Weipeng Xu
# @Date:   2019-02-27 21:57:58
# @Last Modified by:   Weipeng Xu
# @Last Modified time: 2019-02-27 23:15:54

DATA_POS="test_pos.fas"
DATA_NEG="test_neg.fas"

SEED_PARAMS=
SEED_PARAMS="$SEED_PARAMS --min_support 1.0"
SEED_PARAMS="$SEED_PARAMS --word_length 10"
SEED_PARAMS="$SEED_PARAMS --match_file save_motifs.xml"

SEED_PARAMS="$SEED_PARAMS $DATA_POS $DATA_NEG"

# For time and memory usage, please add the following command:
# /usr/bin/time -v  on Linux
# /time -l          on MacOS
../../src/algorithms/seed_w $SEED_PARAMS > ./run.out
