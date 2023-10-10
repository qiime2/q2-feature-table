
export JOBNAME=$parsl.htex.block-0.1696874964.8535182
set -e
export CORES=$(getconf _NPROCESSORS_ONLN)
[[ "1" == "1" ]] && echo "Found cores : $CORES"
WORKERCOUNT=1
FAILONANY=0
PIDS=""

CMD() {
process_worker_pool.py  --max_workers=7 -a isaiah-hplaptop15dy1xxx,134.114.101.49,127.0.0.1,10.18.35.138 -p 0 -c 1.0 -m None --poll 10 --task_port=54950 --result_port=54286 --logdir=/home/isaiahr/work/src/q2-feature-table/runinfo/000/htex --block_id=0 --hb_period=30  --hb_threshold=120 --cpu-affinity none --available-accelerators  --start-method spawn
}
for COUNT in $(seq 1 1 $WORKERCOUNT); do
    [[ "1" == "1" ]] && echo "Launching worker: $COUNT"
    CMD $COUNT &
    PIDS="$PIDS $!"
done

ALLFAILED=1
ANYFAILED=0
for PID in $PIDS ; do
    wait $PID
    if [ "$?" != "0" ]; then
        ANYFAILED=1
    else
        ALLFAILED=0
    fi
done

[[ "1" == "1" ]] && echo "All workers done"
if [ "$FAILONANY" == "1" ]; then
    exit $ANYFAILED
else
    exit $ALLFAILED
fi
