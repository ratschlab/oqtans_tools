#!/bin/sh

numtest=$1
failed=`find -name "*.test_failed" 2>/dev/null|wc -l`
passed=`find -name "*.test_passed" 2>/dev/null|wc -l`

echo "--------------------------------------------------------"
echo "Test summary ($numtest testcases):"
echo "   Test Failed: $failed"
if [ "$failed" -ne "0" ] 
then
	for f in `ls *.test_failed` 
	do
		test=`echo $f | cut -f 1 -d.`
		index=`echo $test | cut -b 1-3`
		if [ "$index" == "bwa" ]
		then
			test=`echo $test | cut -b 4-`
			echo "      $test (bwt-based index): NOK"
		else
			echo "      $test (array-based index): NOK"
		fi
	done
fi
echo "   Test Passed: $passed" 
if [ "$passed" -ne "0" ] 
then
	for f in `ls *.test_passed` 
	do
		test=`echo $f | cut -f 1 -d.`
		index=`echo $test | cut -b 1-3`
		if [ "$index" == "bwa" ]
		then
			test=`echo $test | cut -b 4-`
			echo "      $test (bwt-based index): OK"
		else
			echo "      $test (array-based index): OK"
		fi
	done
fi

rm -f *.test_failed *.test_passed
