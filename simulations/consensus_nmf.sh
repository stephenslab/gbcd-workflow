cnmf prepare --output-dir ./ --name iter2 -c ./iter2.txt -k 9 10 11 12  --n-iter 100 --seed 1  --numgenes 4000
cnmf factorize --output-dir ./ --name iter2 --worker-index 0
cnmf combine --output-dir ./ --name iter2
cnmf k_selection_plot --output-dir ./ --name iter2
cnmf consensus --output-dir ./ --name iter2 --components 12 --local-density-threshold 0.1 --show-clustering


