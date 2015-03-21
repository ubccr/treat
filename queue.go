package treat

type AlignmentQueue []*Alignment

func (pq AlignmentQueue) Len() int { return len(pq) }

func (pq AlignmentQueue) Less(i, j int) bool {
    return pq[i].ReadCount < pq[j].ReadCount
}

func (pq AlignmentQueue) Swap(i, j int) {
    pq[i], pq[j] = pq[j], pq[i]
}

func (pq *AlignmentQueue) Push(x interface{}) {
    item := x.(*Alignment)
    *pq = append(*pq, item)
}

func (pq *AlignmentQueue) Pop() interface{} {
    old := *pq
    n := len(old)
    item := old[n-1]
    *pq = old[0 : n-1]
    return item
}
