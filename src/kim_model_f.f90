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


module kim_model_f_module
  implicit none
  private

  public &
    create, &
    destroy, &
    get_influence_distance, &
    get_cutoffs, &
    get_argument_attribute, &
    get_call_back_attribute, &
    get_units, &
    set_data_int, &
    set_data_double, &
    set_call_back, &
    compute, &
    clear_pointers_and_reinitialize_model, &
    get_species_support_and_code, &
    get_num_params, &
    get_parameter_data_type_and_description, &
    get_parameter_int_extent_and_pointer, &
    get_parameter_double_extent_and_pointer, &
    set_sim_buffer, &
    get_sim_buffer, &
    model_string

  interface
    integer(c_int) function create(numbering, requested_length_unit, &
      requested_energy_unit, requested_charge_unit, &
      requested_temperature_unit, requested_time_unit, model_name, &
      requested_units_accepted, model) bind(c, name="KIM_Model_create")
      use, intrinsic :: iso_c_binding
      use kim_numbering_module, only : kim_numbering_type
      use kim_unit_system_module, only : kim_length_unit_type, &
        kim_energy_unit_type, kim_charge_unit_type, kim_temperature_unit_type, &
        kim_time_unit_type
      implicit none
      type(kim_numbering_type), intent(in), value :: numbering
      type(kim_length_unit_type), intent(in), value :: requested_length_unit
      type(kim_energy_unit_type), intent(in), value :: requested_energy_unit
      type(kim_charge_unit_type), intent(in), value :: requested_charge_unit
      type(kim_temperature_unit_type), intent(in), value :: &
        requested_temperature_unit
      type(kim_time_unit_type), intent(in), value :: requested_time_unit
      character(c_char), intent(in) :: model_name(*)
      integer(c_int), intent(out) :: requested_units_accepted
      type(c_ptr), intent(out) :: model
    end function create

    subroutine destroy(model) bind(c, name="KIM_Model_destroy")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: model
    end subroutine destroy

    subroutine get_influence_distance(model, influence_distance) &
      bind(c, name="KIM_Model_get_influence_distance")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(in) :: model
      real(c_double), intent(out) :: influence_distance
    end subroutine get_influence_distance

    subroutine get_cutoffs(model, number_of_cutoffs, cutoffs_ptr) &
      bind(c, name="KIM_Model_get_cutoffs")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(in) :: model
      integer(c_int), intent(out) :: number_of_cutoffs
      type(c_ptr), intent(out) :: cutoffs_ptr
    end subroutine get_cutoffs

    integer(c_int) function get_argument_attribute(model, argument_name, &
      attribute) bind(c, name="KIM_Model_get_argument_attribute")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      use kim_argument_name_module, only : kim_argument_name_type
      use kim_attribute_module, only : kim_attribute_type
      implicit none
      type(kim_model_type), intent(in) :: model
      type(kim_argument_name_type), intent(in), value :: argument_name
      type(kim_attribute_type), intent(out) :: attribute
    end function get_argument_attribute

    integer(c_int) function get_call_back_attribute(model, call_back_name, &
      attribute) bind(c, name="KIM_Model_get_call_back_attribute")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      use kim_call_back_name_module, only : kim_call_back_name_type
      use kim_attribute_module, only : kim_attribute_type
      implicit none
      type(kim_model_type), intent(in) :: model
      type(kim_call_back_name_type), intent(in), value :: call_back_name
      type(kim_attribute_type), intent(out) :: attribute
    end function get_call_back_attribute

    subroutine get_units(model, length_unit, energy_unit, charge_unit, &
      temperature_unit, time_unit) bind(c, name="KIM_Model_get_units")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      use kim_unit_system_module, only : kim_length_unit_type, &
        kim_energy_unit_type, kim_charge_unit_type, kim_temperature_unit_type, &
        kim_time_unit_type
      type(kim_model_type), intent(in) :: model
      type(kim_length_unit_type), intent(out) :: length_unit
      type(kim_energy_unit_type), intent(out) :: energy_unit
      type(kim_charge_unit_type), intent(out) :: charge_unit
      type(kim_temperature_unit_type), intent(out) :: temperature_unit
      type(kim_time_unit_type), intent(out) :: time_unit
    end subroutine get_units

    integer(c_int) function set_data_int(model, argument_name, ptr) &
      bind(c, name="KIM_Model_set_data_int")
      use, intrinsic :: iso_c_binding
      use kim_argument_name_module, only : kim_argument_name_type
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(inout) :: model
      type(kim_argument_name_type), intent(in), value :: argument_name
      type(c_ptr), intent(in), value :: ptr
    end function set_data_int

    integer(c_int) function set_data_double(model, argument_name, ptr) &
      bind(c, name="KIM_Model_set_data_double")
      use, intrinsic :: iso_c_binding
      use kim_argument_name_module, only : kim_argument_name_type
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(inout) :: model
      type(kim_argument_name_type), intent(in), value :: argument_name
      type(c_ptr), intent(in), value :: ptr
    end function set_data_double

    subroutine set_call_back(model, call_back_name, language_name, &
      fptr, data_object) &
      bind(c, name="KIM_Model_set_call_back")
      use, intrinsic :: iso_c_binding
      use kim_language_name_module, only : kim_language_name_type
      use kim_call_back_name_module, only : kim_call_back_name_type
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(inout) :: model
      type(kim_language_name_type), intent(in), value :: language_name
      type(kim_call_back_name_type), intent(in), value :: call_back_name
      type(c_funptr), intent(in), value :: fptr
      type(c_ptr), intent(in), value :: data_object
    end subroutine set_call_back

    integer(c_int) function compute(model) bind(c, name="KIM_Model_compute")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(in) :: model
    end function compute

    integer(c_int) function clear_pointers_and_reinitialize_model(&
      model) &
      bind(c, name="KIM_Model_clear_pointers_and_reinitialize_model")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(inout) :: model
    end function clear_pointers_and_reinitialize_model

    integer(c_int) function get_species_support_and_code(model, species_name, &
      species_is_supported, code) &
      bind(c, name="KIM_Model_get_species_support_and_code")
      use, intrinsic :: iso_c_binding
      use kim_species_name_module, only : kim_species_name_type
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(in) :: model
      type(kim_species_name_type), intent(in), value :: species_name
      integer(c_int), intent(out) :: species_is_supported
      integer(c_int), intent(out) :: code
    end function get_species_support_and_code

    subroutine get_num_params(model, number_of_parameters) &
      bind(c, name="KIM_Model_get_num_params")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(in) :: model
      integer(c_int), intent(out) :: number_of_parameters
    end subroutine get_num_params

    integer(c_int) function get_parameter_data_type_and_description(model, &
      index, data_type, description) &
      bind(c, name="KIM_Model_get_parameter_data_type_and_description")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      use kim_data_type_module, only : kim_data_type_type
      implicit none
      type(kim_model_type), intent(in) :: model
      integer(c_int), intent(in), value :: index
      type(kim_data_type_type), intent(out) :: data_type
      type(c_ptr), intent(out) :: description
    end function get_parameter_data_type_and_description

    integer(c_int) function get_parameter_int_extent_and_pointer(model, &
      index, extent, ptr) &
      bind(c, name="KIM_Model_get_parameter_int_extent_and_pointer")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(in) :: model
      integer(c_int), intent(in), value :: index
      integer(c_int), intent(out) :: extent
      type(c_ptr), intent(out) :: ptr
    end function get_parameter_int_extent_and_pointer

    integer(c_int) function get_parameter_double_extent_and_pointer(model, &
      index, extent, ptr) &
      bind(c, name="KIM_Model_get_parameter_double_extent_and_pointer")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(in) :: model
      integer(c_int), intent(in), value :: index
      integer(c_int), intent(out) :: extent
      type(c_ptr), intent(out) :: ptr
    end function get_parameter_double_extent_and_pointer

    subroutine set_sim_buffer(model, ptr) &
      bind(c, name="KIM_Model_set_sim_buffer")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(inout) :: model
      type(c_ptr), intent(in), value :: ptr
    end subroutine set_sim_buffer

    subroutine get_sim_buffer(model, ptr) &
      bind(c, name="KIM_Model_get_sim_buffer")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(in) :: model
      type(c_ptr), intent(out) :: ptr
    end subroutine get_sim_buffer

    type(c_ptr) function model_string(model) &
      bind(c, name="KIM_Model_string")
      use, intrinsic :: iso_c_binding
      use kim_model_module, only : kim_model_type
      implicit none
      type(kim_model_type), intent(in) :: model
    end function model_string
  end interface
end module kim_model_f_module


! free functions to implement kim_model_module

subroutine kim_model_create(numbering, requested_length_unit, &
  requested_energy_unit, requested_charge_unit, &
  requested_temperature_unit, requested_time_unit, model_name, &
  requested_units_accepted, model, ierr)
  use, intrinsic :: iso_c_binding
  use kim_numbering_module, only : kim_numbering_type
  use kim_unit_system_module, only : kim_length_unit_type, &
    kim_energy_unit_type, kim_charge_unit_type, kim_temperature_unit_type, &
    kim_time_unit_type
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : create
  implicit none
  type(kim_numbering_type), intent(in), value :: numbering
  type(kim_length_unit_type), intent(in), value :: requested_length_unit
  type(kim_energy_unit_type), intent(in), value :: requested_energy_unit
  type(kim_charge_unit_type), intent(in), value :: requested_charge_unit
  type(kim_temperature_unit_type), intent(in), value :: &
    requested_temperature_unit
  type(kim_time_unit_type), intent(in), value :: requested_time_unit
  character(len=*), intent(in) :: model_name
  integer(c_int), intent(out) :: requested_units_accepted
  type(kim_model_type), intent(out), pointer :: model
  integer(c_int), intent(out) :: ierr

  type(c_ptr) :: pmodel

  ierr = create(numbering, requested_length_unit, requested_energy_unit, &
    requested_charge_unit, requested_temperature_unit, requested_time_unit, &
    trim(model_name)//c_null_char, requested_units_accepted, pmodel)
  call c_f_pointer(pmodel, model)
end subroutine kim_model_create

subroutine kim_model_destroy(model)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : destroy
  implicit none
  type(kim_model_type), intent(inout), pointer :: model

  type(c_ptr) :: pmodel
  pmodel = c_loc(model)
  call destroy(pmodel)
  nullify(model)
end subroutine kim_model_destroy

subroutine kim_model_get_influence_distance(model, influence_distance)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : get_influence_distance
  implicit none
  type(kim_model_type), intent(in) :: model
  real(c_double), intent(out) :: influence_distance

  call get_influence_distance(model, influence_distance)
end subroutine kim_model_get_influence_distance

subroutine kim_model_get_cutoffs(model, number_of_cutoffs, cutoffs)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : get_cutoffs
  implicit none
  type(kim_model_type), intent(in) :: model
  integer(c_int), intent(out) :: number_of_cutoffs
  real(c_double), intent(out), pointer :: cutoffs(:)

  type(c_ptr) cutoffs_ptr

  call get_cutoffs(model, number_of_cutoffs, cutoffs_ptr)
  call c_f_pointer(cutoffs_ptr, cutoffs, [number_of_cutoffs])
end subroutine kim_model_get_cutoffs

subroutine kim_model_get_argument_attribute(model, argument_name, attribute, &
  ierr)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_argument_name_module, only : kim_argument_name_type
  use kim_attribute_module, only : kim_attribute_type
  use kim_model_f_module, only : get_argument_attribute
  implicit none
  type(kim_model_type), intent(in) :: model
  type(kim_argument_name_type), intent(in), value :: argument_name
  type(kim_attribute_type), intent(out) :: attribute
  integer(c_int), intent(out) :: ierr

  ierr = get_argument_attribute(model, argument_name, attribute)
end subroutine kim_model_get_argument_attribute

subroutine kim_model_get_call_back_attribute(model, call_back_name, attribute, &
  ierr)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_call_back_name_module, only : kim_call_back_name_type
  use kim_attribute_module, only : kim_attribute_type
  use kim_model_f_module, only : get_call_back_attribute
  implicit none
  type(kim_model_type), intent(in) :: model
  type(kim_call_back_name_type), intent(in), value :: call_back_name
  type(kim_attribute_type), intent(out) :: attribute
  integer(c_int), intent(out) :: ierr

  ierr = get_call_back_attribute(model, call_back_name, attribute)
end subroutine kim_model_get_call_back_attribute

subroutine kim_model_get_units(model, length_unit, energy_unit, &
  charge_unit, temperature_unit, time_unit)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_unit_system_module, only : kim_length_unit_type, &
    kim_energy_unit_type, kim_charge_unit_type, kim_temperature_unit_type, &
    kim_time_unit_type
  use kim_model_f_module, only : get_units
  type(kim_model_type), intent(in) :: model
  type(kim_length_unit_type), intent(out) :: length_unit
  type(kim_energy_unit_type), intent(out) :: energy_unit
  type(kim_charge_unit_type), intent(out) :: charge_unit
  type(kim_temperature_unit_type), intent(out) :: temperature_unit
  type(kim_time_unit_type), intent(out) :: time_unit

  call get_units(model, length_unit, energy_unit, charge_unit, &
    temperature_unit, time_unit)
end subroutine kim_model_get_units

subroutine kim_model_set_data_int0(model, argument_name, int0, ierr)
  use, intrinsic :: iso_c_binding
  use kim_argument_name_module, only : kim_argument_name_type
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : set_data_int
  implicit none
  type(kim_model_type), intent(inout) :: model
  type(kim_argument_name_type), intent(in), value :: argument_name
  integer(c_int), intent(in), target :: int0
  integer(c_int), intent(out) :: ierr

  ierr = set_data_int(model, argument_name, c_loc(int0))
end subroutine kim_model_set_data_int0

subroutine kim_model_set_data_int1(model, argument_name, int1, ierr)
  use, intrinsic :: iso_c_binding
  use kim_argument_name_module, only : kim_argument_name_type
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : set_data_int
  implicit none
  type(kim_model_type), intent(inout) :: model
  type(kim_argument_name_type), intent(in), value :: argument_name
  integer(c_int), intent(in), target :: int1(:)
  integer(c_int), intent(out) :: ierr

  call set_data(model, argument_name, size(int1,1,c_int), int1, ierr)
  return

contains
  subroutine set_data(model, argument_name, extent1, int1, ierr)
    use, intrinsic :: iso_c_binding
    use kim_argument_name_module, only : kim_argument_name_type
    use kim_model_module, only : kim_model_type
    implicit none
    type(kim_model_type), intent(inout) :: model
    type(kim_argument_name_type), intent(in), value :: argument_name
    integer(c_int), intent(in) :: extent1
    integer(c_int), intent(in), target :: int1(extent1)
    integer(c_int), intent(out) :: ierr

    ierr = set_data_int(model, argument_name, c_loc(int1))
  end subroutine set_data
end subroutine kim_model_set_data_int1

subroutine kim_model_set_data_int2(model, argument_name, int2, ierr)
  use, intrinsic :: iso_c_binding
  use kim_argument_name_module, only : kim_argument_name_type
  use kim_model_module, only : kim_model_type
  implicit none
  type(kim_model_type), intent(inout) :: model
  type(kim_argument_name_type), intent(in), value :: argument_name
  integer(c_int), intent(in), target :: int2(:,:)
  integer(c_int), intent(out) :: ierr

  call set_data(model, argument_name, size(int2, 1, c_int),&
    size(int2, 2, c_int), int2, ierr)
  return

contains
  subroutine set_data(model, argument_name, extent1, extent2, int2, ierr)
    use, intrinsic :: iso_c_binding
    use kim_argument_name_module, only : kim_argument_name_type
    use kim_model_module, only : kim_model_type
    use kim_model_f_module, only : set_data_int
    implicit none
    type(kim_model_type), intent(inout) :: model
    type(kim_argument_name_type), intent(in), value :: argument_name
    integer(c_int), intent(in) :: extent1
    integer(c_int), intent(in) :: extent2
    integer(c_int), intent(in), target :: int2(extent1,extent2)
    integer(c_int), intent(out) :: ierr

    ierr = set_data_int(model, argument_name, c_loc(int2))
  end subroutine set_data
end subroutine kim_model_set_data_int2

subroutine kim_model_set_data_double0(model, argument_name, double0, ierr)
  use, intrinsic :: iso_c_binding
  use kim_argument_name_module, only : kim_argument_name_type
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : set_data_double
  implicit none
  type(kim_model_type), intent(inout) :: model
  type(kim_argument_name_type), intent(in), value :: argument_name
  real(c_double), intent(in), target :: double0
  integer(c_int), intent(out) :: ierr

  ierr = set_data_double(model, argument_name, c_loc(double0))
end subroutine kim_model_set_data_double0

subroutine kim_model_set_data_double1(model, argument_name, double1, ierr)
  use, intrinsic :: iso_c_binding
  use kim_argument_name_module, only : kim_argument_name_type
  use kim_model_module, only : kim_model_type
  implicit none
  type(kim_model_type), intent(inout) :: model
  type(kim_argument_name_type), intent(in), value :: argument_name
  real(c_double), intent(in), target :: double1(:)
  integer(c_int), intent(out) :: ierr

  call set_data(model, argument_name, size(double1, 1, c_int), double1, ierr)
  return

contains
  subroutine set_data(model, argument_name, extent1, double1, ierr)
    use, intrinsic :: iso_c_binding
    use kim_argument_name_module, only : kim_argument_name_type
    use kim_model_module, only : kim_model_type
    use kim_model_f_module, only : set_data_double
    implicit none
    type(kim_model_type), intent(inout) :: model
    type(kim_argument_name_type), intent(in), value :: argument_name
    integer(c_int), intent(in) :: extent1
    real(c_double), intent(in), target :: double1(extent1)
    integer(c_int), intent(out) :: ierr

    ierr = set_data_double(model, argument_name, c_loc(double1))
  end subroutine set_data
end subroutine kim_model_set_data_double1

subroutine kim_model_set_data_double2(model, argument_name, double2, ierr)
  use, intrinsic :: iso_c_binding
  use kim_argument_name_module, only : kim_argument_name_type
  use kim_model_module, only : kim_model_type
  implicit none
  type(kim_model_type), intent(inout) :: model
  type(kim_argument_name_type), intent(in), value :: argument_name
  real(c_double), intent(in), target :: double2(:,:)
  integer(c_int), intent(out) :: ierr

  call set_data(model, argument_name, size(double2, 1, c_int), &
    size(double2, 2, c_int), double2, ierr)
  return

contains
  subroutine set_data(model, argument_name, extent1, extent2, double2, ierr)
    use, intrinsic :: iso_c_binding
    use kim_argument_name_module, only : kim_argument_name_type
    use kim_model_module, only : kim_model_type
    use kim_model_f_module, only : set_data_double
    implicit none
    type(kim_model_type), intent(inout) :: model
    type(kim_argument_name_type), intent(in), value :: argument_name
    integer(c_int), intent(in) :: extent1
    integer(c_int), intent(in) :: extent2
    real(c_double), intent(in), target :: double2(extent1,extent2)
    integer(c_int), intent(out) :: ierr

    ierr = set_data_double(model, argument_name, c_loc(double2))
  end subroutine set_data
end subroutine kim_model_set_data_double2

subroutine kim_model_set_call_back(model, call_back_name, language_name, &
  fptr, data_object)
  use, intrinsic :: iso_c_binding
  use kim_call_back_name_module, only : kim_call_back_name_type
  use kim_language_name_module, only : kim_language_name_type
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : set_call_back
  implicit none
  type(kim_model_type), intent(inout) :: model
  type(kim_call_back_name_type), intent(in), value :: call_back_name
  type(kim_language_name_type), intent(in), value :: language_name
  type(c_funptr), intent(in), value :: fptr
  type(c_ptr), intent(in), value :: data_object

  call set_call_back(model, call_back_name, language_name, fptr, data_object)
end subroutine kim_model_set_call_back

subroutine kim_model_compute(model, ierr)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : compute
  implicit none
  type(kim_model_type), intent(in) :: model
  integer(c_int), intent(out) :: ierr

  ierr = compute(model)
end subroutine kim_model_compute

subroutine kim_model_clear_pointers_and_reinitialize_model(model, ierr)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : clear_pointers_and_reinitialize_model
  implicit none
  type(kim_model_type), intent(inout) :: model
  integer(c_int), intent(out) :: ierr

  ierr = clear_pointers_and_reinitialize_model(model)
end subroutine kim_model_clear_pointers_and_reinitialize_model

subroutine kim_model_get_species_support_and_code(model, species_name, &
  species_is_supported, code, ierr)
  use, intrinsic :: iso_c_binding
  use kim_species_name_module, only : kim_species_name_type
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : get_species_support_and_code
  implicit none
  type(kim_model_type), intent(in) :: model
  type(kim_species_name_type), intent(in), value :: species_name
  integer(c_int), intent(out) :: species_is_supported
  integer(c_int), intent(out) :: code
  integer(c_int), intent(out) :: ierr

  ierr = get_species_support_and_code(model, species_name, &
    species_is_supported, code)
end subroutine kim_model_get_species_support_and_code

subroutine kim_model_get_num_params(model, number_of_parameters)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : get_num_params
  implicit none
  type(kim_model_type), intent(in) :: model
  integer(c_int), intent(out) :: number_of_parameters

  call get_num_params(model, number_of_parameters)
end subroutine kim_model_get_num_params

subroutine kim_model_get_parameter_data_type_and_description(model, index, &
  data_type, description, ierr)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_data_type_module, only : kim_data_type_type
  use kim_model_f_module, only : get_parameter_data_type_and_description
  implicit none
  type(kim_model_type), intent(in) :: model
  integer(c_int), intent(in), value :: index
  type(kim_data_type_type), intent(out) :: data_type
  character(len=*), intent(out) :: description
  integer(c_int), intent(out) :: ierr

  type(c_ptr) :: p
  character(len=len(description)), pointer :: fp
  integer(c_int) :: null_index

  ierr = get_parameter_data_type_and_description(model, index-1, data_type, p)
  call c_f_pointer(p, fp)
  null_index = scan(fp, char(0))-1
  description = fp(1:null_index)
end subroutine kim_model_get_parameter_data_type_and_description

subroutine kim_model_get_parameter_int_extent_and_pointer(model, index, &
  extent, int1, ierr)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : get_parameter_int_extent_and_pointer
  implicit none
  type(kim_model_type), intent(in) :: model
  integer(c_int), intent(in), value :: index
  integer(c_int), intent(out) :: extent
  integer(c_int), intent(out), pointer :: int1(:)
  integer(c_int), intent(out) :: ierr

  type(c_ptr) p

  ierr = get_parameter_int_extent_and_pointer(model, index-1, extent, p)
  call c_f_pointer(p, int1, [extent])
end subroutine kim_model_get_parameter_int_extent_and_pointer

subroutine kim_model_get_parameter_double_extent_and_pointer(model, index, &
  extent, double1, ierr)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : get_parameter_double_extent_and_pointer
  implicit none
  type(kim_model_type), intent(in) :: model
  integer(c_int), intent(in), value :: index
  integer(c_int), intent(out) :: extent
  real(c_double), intent(out), pointer :: double1(:)
  integer(c_int), intent(out) :: ierr

  type(c_ptr) p

  ierr = get_parameter_double_extent_and_pointer(model, index-1, extent, p)
  call c_f_pointer(p, double1, [extent])
end subroutine kim_model_get_parameter_double_extent_and_pointer

subroutine kim_model_set_sim_buffer(model, ptr)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : set_sim_buffer
  implicit none
  type(kim_model_type), intent(inout) :: model
  type(c_ptr), intent(in), value :: ptr

  call set_sim_buffer(model, ptr)
end subroutine kim_model_set_sim_buffer

subroutine kim_model_get_sim_buffer(model, ptr)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : get_sim_buffer
  implicit none
  type(kim_model_type), intent(in) :: model
  type(c_ptr), intent(out) :: ptr

  call get_sim_buffer(model, ptr)
end subroutine kim_model_get_sim_buffer

subroutine kim_model_string(model, string)
  use, intrinsic :: iso_c_binding
  use kim_model_module, only : kim_model_type
  use kim_model_f_module, only : model_string
  implicit none
  type(kim_model_type), intent(in) :: model
  character(len=*), intent(out) :: string

  type(c_ptr) :: p
  character(len=len(string)), pointer :: fp
  integer(c_int) :: null_index

  p = model_string(model)
  call c_f_pointer(p, fp)
  null_index = scan(fp, char(0))-1
  string = fp(1:null_index)
end subroutine kim_model_string
