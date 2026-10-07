#!/usr/bin/env bash
cd "$(dirname "${BASH_SOURCE[0]}")"
target="plots-for-paper"

# The tissue.png images probably won't be used in the paper,
# but we collect them here anyway because they are a useful reference
# when inspecting the transcripts.png images.

mkdir -p "$target/BMKMANU S1000"
cp "plots/BMKMANU S1000/exp74/analysis/tissue.png" "$target/BMKMANU S1000/"
cp "plots/BMKMANU S1000/exp74/analysis/transcripts.png" "$target/BMKMANU S1000/"

mkdir -p "$target/BMKMANU S3000"
cp "plots/BMKMANU S3000/exp99/analysis/tissue.png" "$target/BMKMANU S3000/"
cp "plots/BMKMANU S3000/exp99/analysis/transcripts.png" "$target/BMKMANU S3000/"

mkdir -p "$target/Nova-ST"
cp "plots/Nova-ST/Exp65-sampleSPC014/analysis/tissue.png" "$target/Nova-ST/"
cp "plots/Nova-ST/Exp65-sampleSPC014/analysis/transcripts.png" "$target/Nova-ST/"

mkdir -p "$target/Visium HD (FFPE)"
cp "plots/Visium HD (FFPE)/SCA003/analysis/tissue.png" "$target/Visium HD (FFPE)/"
cp "plots/Visium HD (FFPE)/SCA003/analysis/transcripts.png" "$target/Visium HD (FFPE)/"

mkdir -p "$target/Visium HD (FF)"
cp "plots/Visium HD (FF)/SCA019/analysis/tissue.png" "$target/Visium HD (FF)/"
cp "plots/Visium HD (FF)/SCA019/analysis/transcripts.png" "$target/Visium HD (FF)/"

mkdir -p "$target/Visium HD 3'"
cp "plots/Visium HD 3'/SCA024/analysis/tissue.png" "$target/Visium HD 3'/"
cp "plots/Visium HD 3'/SCA024/analysis/transcripts.png" "$target/Visium HD 3'/"

mkdir -p "$target/Stereo-Seq"
cp "plots/Stereo-Seq/spc014/analysis/tissue.png" "$target/Stereo-Seq/"
cp "plots/Stereo-Seq/spc014/analysis/transcripts.png" "$target/Stereo-Seq/"

cp "plots/combined/cdf_distance.png" "$target/"
cp "plots/combined/binwidth_1um/distance_kde_um_peaksaligned.png" "$target/"
