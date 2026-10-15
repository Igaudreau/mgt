#!/bin/bash
#SBATCH --mem=4G
#SBATCH --time=0-03:00
#SBATCH --cpus-per-task=1
#SBATCH --account=def-bmaxwell
#SBATCH --output=movieout

# First file to start from
START=1

for i in rho te rgrad sch; do
  rm $i.mp4

  if hash ffmpeg 2>/dev/null; then
    ffmpeg -start_number $START -i ./png/${i}%04d.png -c:v libx264 -crf 23 -profile:v baseline -level 3.0 -pix_fmt yuv420p -c:a aac -ac 2 -b:a 128k -movflags faststart -s 1920x1080 $i.mp4

 # elif hash avconv 2>/dev/null; then
 # avconv -f image2 -start_number 00001 -i ./png/${i}%04d.png -r 10 -s 640x480 $i.avi

  else
    echo "movie compiler (ffmpeg) not found!"
  fi

done

if hash ffmpeg 2>/dev/null; then
  rm collage.avi
  ffmpeg -i sch.mp4 -i te.mp4 -i rho.mp4 -i rgrad.mp4 -vcodec libx264 -filter_complex "[0:v][1:v]hstack=inputs=2[top]; [2:v][3:v]hstack=inputs=2[bottom]; [top][bottom]vstack=inputs=2[v]; [v]scale=-1:1500[out]" -map "[out]" -s 1920x1080 collage.avi
else
  echo "movie compiler (ffmpeg) not found!"
fi


