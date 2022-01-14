# GEB-2021-0553

Here, we use [`Docker`]() technology for creating a virtual system to reporoduce the results in our analysis. Docker containers are very common for creating micro-services...


# Installation
Docker can be installed in different platforms. You can find the instruction for installing docker on different operating systems visit [docker website](https://docs.docker.com/get-docker/).

# Creating the image
Follow the inststruction below to create the docker image with RSudio installed and all the R and system packages required for running our analysis.

The following command are terminal commands.

## Clone this repository
You need to have [`git`]() installed on your system.

```bash
git clone https://github.com/rvalavi/GEB-2021-0553.git
```

## Build the docker image
In Linux systmes you might need to use `sudo` command before docker.

```bash
cd GEB-2021-0553

docker build -t rvalavi:4.0 .
```

Wait until the build is complete.

## Run the docker container


```bash
docker run --name rstudio -p 8787:8787 -e PASSWORD=123 -v /home/rvalavi:project -d rvalavi:4.0

```

## Run RStudio server
Open an internet browser and go to `localhost:8787` to open RStudio. Use the password you specified in the previous step (here is: 123) to open RStudio.
