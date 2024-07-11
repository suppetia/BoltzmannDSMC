module m_quadtree
  use m_types, only: fp, pp, i4, dp
  implicit none

  type :: QuadTreeNode
    integer(pp) :: value = 0

    !> it should be possible to calculate this on the fly
    !> but I'm not capable of wrapping my head around this
    real(dp) :: x, y
    real(dp) :: width, height

    type(QuadTreeNode), allocatable, dimension(:) :: children

  end type QuadTreeNode

  type :: QuadTree
    integer :: maxDepth

    type(QuadTreeNode) :: root = QuadTreeNode(x=1, y=1, width=0, height=0)

  end type QuadTree

contains

  !> add children to a node
  subroutine splitNode(node)
    implicit none
    type(QuadTreeNode), intent(inout) :: node

    integer(i4) :: newWidth, newHeight

    allocate(node%children(4)) !> all children now have the value 0
    !> handle this accordingly

    newWidth = node%width / 2
    newHeight = node%height / 2
    
    node%children(1) = QuadTreeNode(x=node%x, y=node%y, width=newWidth, height=newHeight)
    node%children(2) = QuadTreeNode(x=node%x+newWidth, y=node%y, width=newWidth, height=newHeight)
    node%children(3) = QuadTreeNode(x=node%x, y=node%y+newHeight, width=newWidth, height=newHeight)
    node%children(4) = QuadTreeNode(x=node%x+newWidth, y=node%y+newHeight, width=newWidth, height=newHeight)
    
    !> mark that the node is no leave
    node%value = -1

  end subroutine splitNode

  !> merge the children making node a leaf
  subroutine mergeChildNodes(node)
    implicit none
    type(QuadTreeNode), intent(inout) :: node

    !> TODO: handle the value of node

    deallocate(node%children)
  end subroutine mergeChildNodes

  !> delete a QuadTree by freeing all memory
  recursive subroutine deleteSubTree(node)
    implicit none
    type(QuadTreeNode), intent(inout) :: node
    integer :: i
    
    if (allocated(node%children)) then
      do i = 1,4
        call deleteSubTree(node%children(i))
      end do
      deallocate(node%children)
    end if
  end subroutine deleteSubTree
    


end module m_quadtree


