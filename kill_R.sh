ps -u ashtonsb | grep R | awk '{print $1}' | xargs kill -KILL
