for f in $(docker exec als bash -c "ls Rplots*")
do
	docker cp als:/app/$f ./$f
done
