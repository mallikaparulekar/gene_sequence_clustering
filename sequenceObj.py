class Sequence(object):
    trueCluster = -1
    value = []
    currentCluster = -1


    def make_sequence(trueCluster, value):
        sequence = Sequence()
        sequence.trueCluster = trueCluster
        sequence.value = value
        return sequence
