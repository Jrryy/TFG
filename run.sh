docker rm als

docker build -t als .

docker run -v useful_files:/app --name als -it als bash
