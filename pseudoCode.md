## Pseudocode of a greedy search for a necklace cover

**ALG:** *Find_One_Necklace_Cover*

Preliminary: Each node has an unique identifier integer in [1,deg(V)]
~~~
INPUT: Adjacency List(adjList) = <node, <neighbors>>
OUTPUT: One necklace cover(not necessarily optimal)

MAIN():
    <paths> = FindPaths(adjList)
    PRINT <paths>

GreedyPath(start, adjList):
    DEFINE vector <path>
    current = start
    ADD "current" to <path>

    WHILE(TRUE):
        tmp = false
        FOR neighbor in adjList:
            IF(neighbor is not visited):
                ADD "neighbor" to <path>
                current = neighbor
                tmp = true
                BREAK
            ENDIF
        ENDFOR
        IF(!tmp) BREAK ENDIF
    ENDWHILE

    RETURN <path>

FindPaths(adjList):
    DEFINE vector <paths>

    FOR node in adjList:
        IF(node is not visited):
            <path> = GreedyPath(node, adjList)
            ADD <path> to <paths>
        ENDIF
    ENDFOR

    RETURN <paths>
~~~