#!/bin/bash


rsync \
   --verbose \
   --progress \
   --archive \
   --delete \
   --delete-excluded \
   --filter="+ /README.md" \
   --filter="+ /.gitignore" \
   --filter="+ /.Rprofile" \
   --filter="+ /renv" \
   --filter="+ /renv/activate.R" \
   --filter="+ /renv.lock" \
   --filter="+ /submission" \
   --filter="+ /submission/**" \
   --filter="+ /src" \
   --filter="+ /src/**" \
   --filter="+ /output" \
   --filter="+ /output/**" \
   --filter="- /**" \
   ahlmanne@seneca:~/projects/perturbation_prediction-benchmark/ benchmark