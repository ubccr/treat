package treat

import (
    "testing"
    "container/heap"
)

func TestQueue(t *testing.T) {

    counts := []uint64{20, 10, 40, 50, 90, 80, 70, 100}
    pq := make(AlignmentQueue, 0)
    heap.Init(&pq)

    for _, v := range(counts) {
        heap.Push(&pq, &Alignment{ReadCount: v})
        if pq.Len() > 4 {
            heap.Pop(&pq)
        }
    }

    if pq.Len() != 4 {
        t.Errorf("Wrong length for priority queue")
    }

    good := []uint64{100, 80, 90, 70}
    for i := pq.Len()-1; i >= 0; i-- {
        if pq[i].ReadCount != good[3-i] {
            t.Errorf("Wrong value for item %d != %d", pq[i].ReadCount, good[3-i])
        }
    }
}
