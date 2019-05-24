docker rm als

docker build -t als .

docker run -d -v useful_files:/app --name als -it als
