!
! CDDL HEADER START
!
! The contents of this file are subject to the terms of the Common Development
! and Distribution License Version 1.0 (the "License").
!
! You can obtain a copy of the license at
! http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
! specific language governing permissions and limitations under the License.
!
! When distributing Covered Code, include this CDDL HEADER in each file and
! include the License file in a prominent location with the name LICENSE.CDDL.
! If applicable, add the following below this CDDL HEADER, with the fields
! enclosed by brackets "[]" replaced with your own identifying information:
!
! Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
!
! CDDL HEADER END
!

!
! Copyright (c) 2016--2017, Regents of the University of Minnesota.
! All rights reserved.
!
! Contributors:
!    Ryan S. Elliott
!

!
! Release: This file is part of the kim-api.git repository.
!


module kim_model_compute_module
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public &
    kim_model_compute_type, &
    kim_model_compute_get_neigh, &
    kim_model_compute_process_dedr, &
    kim_model_compute_process_d2edr2, &
    kim_model_compute_get_data, &
    kim_model_compute_is_call_back_present, &
    kim_model_compute_get_model_buffer, &
    kim_model_compute_log, &
    kim_model_compute_string

  type, bind(c) :: kim_model_compute_type
    private
    type(c_ptr) :: p
  end type kim_model_compute_type

  interface kim_model_compute_get_data
    subroutine kim_model_compute_get_data_int0(model_compute, argument_name, &
      int0, ierr)
      use, intrinsic :: iso_c_binding
      use kim_argument_name_module, only : kim_argument_name_type
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      type(kim_argument_name_type), intent(in), value :: argument_name
      integer(c_int), intent(out), pointer :: int0
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_get_data_int0

    subroutine kim_model_compute_get_data_int1(model_compute, argument_name, &
      extent1, int1, ierr)
      use, intrinsic :: iso_c_binding
      use kim_argument_name_module, only : kim_argument_name_type
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      type(kim_argument_name_type), intent(in), value :: argument_name
      integer(c_int), intent(in), value :: extent1
      integer(c_int), intent(out), pointer :: int1(:)
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_get_data_int1

    subroutine kim_model_compute_get_data_int2(model_compute, argument_name, &
      extent1, extent2, int2, ierr)
      use, intrinsic :: iso_c_binding
      use kim_argument_name_module, only : kim_argument_name_type
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      type(kim_argument_name_type), intent(in), value :: argument_name
      integer(c_int), intent(in), value :: extent1
      integer(c_int), intent(in), value :: extent2
      integer(c_int), intent(out), pointer :: int2(:,:)
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_get_data_int2

    subroutine kim_model_compute_get_data_double0(model_compute, &
      argument_name, double0, ierr)
      use, intrinsic :: iso_c_binding
      use kim_argument_name_module, only : kim_argument_name_type
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      type(kim_argument_name_type), intent(in), value :: argument_name
      real(c_double), intent(out), pointer :: double0
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_get_data_double0

    subroutine kim_model_compute_get_data_double1(model_compute, &
      argument_name, extent1, double1, ierr)
      use, intrinsic :: iso_c_binding
      use kim_argument_name_module, only : kim_argument_name_type
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      type(kim_argument_name_type), intent(in), value :: argument_name
      integer(c_int), intent(in), value :: extent1
      real(c_double), intent(out), pointer :: double1(:)
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_get_data_double1

    subroutine kim_model_compute_get_data_double2(model_compute, &
      argument_name, extent1, extent2, double2, ierr)
      use, intrinsic :: iso_c_binding
      use kim_argument_name_module, only : kim_argument_name_type
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      type(kim_argument_name_type), intent(in), value :: argument_name
      integer(c_int), intent(in), value :: extent1
      integer(c_int), intent(in), value :: extent2
      real(c_double), intent(out), pointer :: double2(:,:)
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_get_data_double2
  end interface kim_model_compute_get_data

  interface
    subroutine kim_model_compute_get_neigh(model_compute, neighbor_list_index, &
      particle_number, number_of_neighbors, neighbors_of_particle, ierr)
      use, intrinsic :: iso_c_binding
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      integer(c_int), intent(in), value :: neighbor_list_index
      integer(c_int), intent(in), value :: particle_number
      integer(c_int), intent(out) :: number_of_neighbors
      integer(c_int), intent(out), pointer :: neighbors_of_particle(:)
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_get_neigh

    subroutine kim_model_compute_process_dedr(model_compute, de, r, dx, i, j, &
      ierr)
      use, intrinsic :: iso_c_binding
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      real(c_double), intent(in), value :: de
      real(c_double), intent(in), value :: r
      type(c_ptr), intent(in) :: dx
      integer(c_int), intent(in), value :: i
      integer(c_int), intent(in), value :: j
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_process_dedr

    subroutine kim_model_compute_process_d2edr2(model_compute, de, r, dx, i, &
      j, ierr)
      use, intrinsic :: iso_c_binding
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      real(c_double), intent(in), value :: de
      type(c_ptr), intent(in), value :: r
      type(c_ptr), intent(in), value :: dx
      type(c_ptr), intent(in), value :: i
      type(c_ptr), intent(in), value :: j
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_process_d2edr2

    subroutine kim_model_compute_is_call_back_present(model_compute, &
      call_back_name, present, ierr)
      use, intrinsic :: iso_c_binding
      use kim_call_back_name_module, only : kim_call_back_name_type
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      type(kim_call_back_name_type), intent(in), value :: call_back_name
      integer(c_int), intent(out) :: present
      integer(c_int), intent(out) :: ierr
    end subroutine kim_model_compute_is_call_back_present

    subroutine kim_model_compute_get_model_buffer(model_compute, ptr)
      use, intrinsic :: iso_c_binding
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      type(c_ptr), intent(out) :: ptr
    end subroutine kim_model_compute_get_model_buffer

    subroutine kim_model_compute_log(model_compute, log_level, message, &
      line_number, file_name)
      use, intrinsic :: iso_c_binding
      use kim_log_level_module, only : kim_log_level_type
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      type(kim_log_level_type), intent(in), value :: log_level
      character(len=*), intent(in) :: message
      integer(c_int), intent(in), value :: line_number
      character(len=*), intent(in) :: file_name
    end subroutine kim_model_compute_log

    subroutine kim_model_compute_string(model_compute, string)
      use, intrinsic :: iso_c_binding
      import kim_model_compute_type
      implicit none
      type(kim_model_compute_type), intent(in) :: model_compute
      character(len=*), intent(out) :: string
    end subroutine kim_model_compute_string
end interface
end module kim_model_compute_module
