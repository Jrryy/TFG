for f in $(docker exec als bash -c "ls Rplot*")
do
	docker cp als:/app/$f ./$f
done
