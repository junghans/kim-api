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


module kim_model_destroy_f_module
  implicit none
  private

  public &
    get_model_buffer, &
    log, &
    model_destroy_string

  interface
    subroutine get_model_buffer(model_destroy, ptr) &
      bind(c, name="KIM_ModelDestroy_get_model_buffer")
      use, intrinsic :: iso_c_binding
      use kim_model_destroy_module, only : kim_model_destroy_type
      implicit none
      type(kim_model_destroy_type), intent(in) :: model_destroy
      type(c_ptr), intent(out) :: ptr
    end subroutine get_model_buffer

    subroutine log(model_destroy, log_level, message, line_number, &
      file_name) bind(c, name="KIM_ModelDestroy_Log")
      use, intrinsic :: iso_c_binding
      use kim_model_destroy_module, only : kim_model_destroy_type
      use kim_log_level_module, only : kim_log_level_type
      implicit none
      type(kim_model_destroy_type), intent(in) :: model_destroy
      type(kim_log_level_type), intent(in), value :: log_level
      character(c_char), intent(in) :: message(*)
      integer(c_int), intent(in), value :: line_number
      character(c_char), intent(in) :: file_name(*)
    end subroutine log

    type(c_ptr) function model_destroy_string(model_destroy) &
      bind(c, name="KIM_ModelDestroy_string")
      use, intrinsic :: iso_c_binding
      use kim_model_destroy_module, only : kim_model_destroy_type
      implicit none
      type(kim_model_destroy_type), intent(in) :: model_destroy
    end function model_destroy_string
  end interface
end module kim_model_destroy_f_module


! free functions to implement kim_model_destroy_module

subroutine kim_model_destroy_get_model_buffer(model_destroy, ptr)
  use, intrinsic :: iso_c_binding
  use kim_model_destroy_module, only : kim_model_destroy_type
  use kim_model_destroy_f_module, only : get_model_buffer
  implicit none
  type(kim_model_destroy_type), intent(in) :: model_destroy
  type(c_ptr), intent(out) :: ptr

  call get_model_buffer(model_destroy, ptr)
end subroutine kim_model_destroy_get_model_buffer

subroutine kim_model_destroy_log(model_destroy, log_level, message, &
  line_number, file_name)
  use, intrinsic :: iso_c_binding
  use kim_model_destroy_module, only : kim_model_destroy_type
  use kim_log_level_module, only : kim_log_level_type
  use kim_model_destroy_f_module, only : log
  implicit none
  type(kim_model_destroy_type), intent(in) :: model_destroy
  type(kim_log_level_type), intent(in), value :: log_level
  character(len=*), intent(in) :: message
  integer(c_int), intent(in), value :: line_number
  character(len=*), intent(in) :: file_name

  call log(model_destroy, log_level, trim(message)//c_null_char, line_number, &
    trim(file_name)//c_null_char)
end subroutine kim_model_destroy_log

subroutine kim_model_destroy_string(model_destroy, string)
  use, intrinsic :: iso_c_binding
  use kim_model_destroy_module, only : kim_model_destroy_type
  use kim_model_destroy_f_module, only : model_destroy_string
  implicit none
  type(kim_model_destroy_type), intent(in) :: model_destroy
  character(len=*), intent(out) :: string

  type(c_ptr) :: p
  character(len=len(string)), pointer :: fp
  integer(c_int) :: null_index

  p = model_destroy_string(model_destroy)
  call c_f_pointer(p, fp)
  null_index = scan(fp, char(0))-1
  string = fp(1:null_index)
end subroutine kim_model_destroy_string
