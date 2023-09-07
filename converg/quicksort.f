      SUBROUTINE quicksort(n,arr)
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
c          Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n
c          is input; arr is replaced on output by its sorted rearrangement.
c          Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
c          auxiliary storage.
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then                      ! Insertion sort when subarray small enough.
           do 12 j=l+1,ir
                a=arr(j)
                do 11 i=j-1,l,-1
                     if(arr(i).le.a)goto 2
                     arr(i+1)=arr(i)
11              enddo
                i=l-1
2               arr(i+1)=a
12         enddo
      
          if(jstack.eq.0)return
          ir=istack(jstack)                ! Pop stack and begin a new round of partitioning.
          l=istack(jstack-1)
          jstack=jstack-2
      else
          k=(l+ir)/2                       ! Choose median of left, center, and right elements
          temp=arr(k)                      ! as partitioning element a. Also rearrange so that
          arr(k)=arr(l+1)                  ! a(l) <= a(l+1) <= a(ir).
          arr(l+1)=temp
          if(arr(l).gt.arr(ir))then
              temp=arr(l)
              arr(l)=arr(ir)
              arr(ir)=temp
          endif
          if(arr(l+1).gt.arr(ir))then
              temp=arr(l+1)
              arr(l+1)=arr(ir)
              arr(ir)=temp
          endif
          if(arr(l).gt.arr(l+1))then
              temp=arr(l)
              arr(l)=arr(l+1)
              arr(l+1)=temp
          endif
          i=l+1                            ! Initialize pointers for partitioning.
          j=ir
          a=arr(l+1)                       ! Partitioning element.
3         continue                         ! Beginning of innermost loop.
              i=i+1                        ! Scan up to find element > a.
          if(arr(i).lt.a)goto 3
4         continue
              j=j-1                        ! Scan down to find element < a.
          if(arr(j).gt.a)goto 4
          if(j.lt.i)goto 5                 ! Pointers crossed. Exit with partitioning complete.
          temp=arr(i)                      ! Exchange elements.
          arr(i)=arr(j)
          arr(j)=temp
          goto 3                           ! End of innermost loop.
5         arr(l+1)=arr(j)                  ! Insert partitioning element.
          arr(j)=a
          jstack=jstack+2
c            Push pointers to larger subarray on stack, process smaller subarray immediately.
          if(jstack.gt.NSTACK) pause "NSTACK too small in sort"
          if(ir-i+1.ge.j-l)then
              istack(jstack)=ir
              istack(jstack-1)=i
              ir=j-1
          else
              istack(jstack)=j-1
              istack(jstack-1)=l
              l=i
          endif
      endif
      goto 1
      END


