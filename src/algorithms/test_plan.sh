# @Author: Weipeng Xu
# @Date:   2019-02-27 21:57:58
# @Last Modified by:   Weipeng Xu
# @Last Modified time: 2019-03-30 11:00:14

str_len=(1000 2000 4000 8000 16000 32000)
str_num=(100 200 400 800)
range=(10 20 40 80 160 320 640)
word_num=(1 2 3 4 6 8)

for i in "${str_len[@]}"
do
DATA_POS="./test_files/str_len/${i}_pos.fas"
DATA_NEG="./test_files/str_len/${i}_neg.fas"
PARAMS="$DATA_POS $DATA_NEG"
SEED_PARAMS=
SEED_PARAMS="$SEED_PARAMS --min_support 0.6"
SEED_PARAMS="$SEED_PARAMS --word_length 10"
SEED_PARAMS="$SEED_PARAMS --match_file ./test_files/str_len/save_motifs.xml"

SEED_PARAMS="$SEED_PARAMS $DATA_POS $DATA_NEG"
/usr/bin/time -o ./test_files/str_len/${i}.time -v ./seed_w $SEED_PARAMS > ./run.out
done

for j in "${range[@]}"
do 
DATA_POS="./test_files/range/${j}_pos.fas"
DATA_NEG="./test_files/range/${j}_neg.fas"
PARAMS="$DATA_POS $DATA_NEG"
SEED_PARAMS=
SEED_PARAMS="$SEED_PARAMS --min_support 0.6"
SEED_PARAMS="$SEED_PARAMS --word_length 10"
SEED_PARAMS="$SEED_PARAMS --match_file ./test_files/range/save_motifs.xml"

SEED_PARAMS="$SEED_PARAMS $DATA_POS $DATA_NEG"
/usr/bin/time -o ./test_files/range/${j}.time -v ./seed_w $SEED_PARAMS > ./run.out
done

for k in "${str_num[@]}"
do 
DATA_POS="./test_files/str_num/${k}_pos.fas"
DATA_NEG="./test_files/str_num/${k}_neg.fas"
PARAMS="$DATA_POS $DATA_NEG"
SEED_PARAMS=
SEED_PARAMS="$SEED_PARAMS --min_support 0.6"
SEED_PARAMS="$SEED_PARAMS --word_length 10"
SEED_PARAMS="$SEED_PARAMS --match_file ./test_files/str_num/save_motifs.xml"

SEED_PARAMS="$SEED_PARAMS $DATA_POS $DATA_NEG"
/usr/bin/time -o ./test_files/str_num/${k}.time -v ./seed_w $SEED_PARAMS > ./run.out
done

for l in "${word_num[@]}"
do 
DATA_POS="./test_files/word_num/${l}_pos.fas"
DATA_NEG="./test_files/word_num/${l}_neg.fas"
PARAMS="$DATA_POS $DATA_NEG"
SEED_PARAMS=
SEED_PARAMS="$SEED_PARAMS --min_support 0.6"
SEED_PARAMS="$SEED_PARAMS --word_length 10"
SEED_PARAMS="$SEED_PARAMS --match_file ./test_files/word_num/save_motifs.xml"

SEED_PARAMS="$SEED_PARAMS $DATA_POS $DATA_NEG"
/usr/bin/time -o ./test_files/word_num/${l}.time -v ./seed_w $SEED_PARAMS > ./run.out
done


DATA_POS="./test_files/limit_test/test_pos.fas"
DATA_NEG="./test_files/limit_test/test_neg.fas"
PARAMS="$DATA_POS $DATA_NEG"
SEED_PARAMS=
SEED_PARAMS="$SEED_PARAMS --min_support 0.6"
SEED_PARAMS="$SEED_PARAMS --word_length 10"
SEED_PARAMS="$SEED_PARAMS --match_file ./test_files/limit_test/save_motifs.xml"

SEED_PARAMS="$SEED_PARAMS $DATA_POS $DATA_NEG"
/usr/bin/time -o ./test_files/limit_test/limit_test.time -v ./seed_w $SEED_PARAMS > ./run.out


# For time and memory usage, please add the following command:
# /usr/bin/time -v  on Linux
# /time -l          on MacOS
# ../../src/algorithms/seed_w $SEED_PARAMS > ./run.out
