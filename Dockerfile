# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install Nextflow
RUN apt-get update && apt-get install -y default-jre curl && \
    curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

# Make port 80 available to the world outside this container
EXPOSE 80

# run the pipeline when the container launches
CMD ["nextflow", "run", "main.nf"]
