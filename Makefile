test:
	R --vanilla -e "devtools::test()" | tee logs/tests.txt
