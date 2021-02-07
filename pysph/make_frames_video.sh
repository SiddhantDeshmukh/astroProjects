cd ./figs/frames-star
ffmpeg -r 30 -f image2 -s 800x600 -i sph-%03d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" evolution.mp4
mv evolution.mp4 ../videos/
cd ../../