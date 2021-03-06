#!/bin/sh
#

#
# CDDL HEADER START
#
# The contents of this file are subject to the terms of the Common Development
# and Distribution License Version 1.0 (the "License").
#
# You can obtain a copy of the license at
# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
# specific language governing permissions and limitations under the License.
#
# When distributing Covered Code, include this CDDL HEADER in each file and
# include the License file in a prominent location with the name LICENSE.CDDL.
# If applicable, add the following below this CDDL HEADER, with the fields
# enclosed by brackets "[]" replaced with your own identifying information:
#
# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2017--2018, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Ryan S. Elliott
#    John S. Spear
#

#
# Release: This file is part of the kim-api.git repository.
#

collections_info=###COLLECTIONS#INFO#UTILITY###
build_config=###BUILD#CONFIG###
major_version=###MAJOR#VERSION###

make_command="make --no-print-directory"


# define usage function
usage () {
  local command=`printf $0 | sed 's|.*/\([^/][^/]*\)/*|\1|'`

  # Follows docopt.org format
  printf "Usage:\n"
  printf "  ${command} list [--with-version]\n"
  printf "  ${command} set-user-models-dir <directory>\n"
  printf "  ${command} set-user-drivers-dir <directory>\n"
  printf "  ${command} install\n"
  printf "          (CWD | environment | user | system [--sudo])\n"
  printf "          (<openkim-item-id>... | <local-item-path>... | OpenKIM)\n"
  printf "  ${command} reinstall [--sudo] (<openkim-item-id> | <local-item-path>)...\n"
  printf "  ${command} remove [--sudo] <item-id>...\n"
  printf "  ${command} remove-all [--sudo]\n"
  printf "\n\n"

  printf "list:\n"
  printf "  List installed kim-api models and model drivers\n"
  printf "\n"
  printf "set-user-models-dir:\n"
  printf "  Rewrite configuration file with provided directory\n"
  printf "\n"
  printf "set-user-drivers-dir:\n"
  printf "  Rewrite configuration file with provided directory\n"
  printf "\n"
  printf "install:\n"
  printf "  Install model and/or model driver from openkim.org or from a local path\n"
  printf "  (Installing to the environment collection places items in the first\n"
  printf "   directory of the list.)\n"
  printf "\n"
  printf "reinstall:\n"
  printf "  Remove and reinstall items.\n"
  printf "\n"
  printf "remove:\n"
  printf "  Remove model or model driver.\n"
  printf "  WARNING: This will remove the entire directory of files associated\n"
  printf "           with the item.\n"
  printf "\n"
  printf "remove-all:\n"
  printf "  Remove all items from all collections.\n"
  printf "  WARNING: This will remove the entire directory of files associated\n"
  printf "           with all items.\n"
}

check_version_compatibility () {
  local version="$1"
  local major=`printf -- "${version}" | sed -e 's/\([^.}]*\).*/\1/'`
  local minor=`printf -- "${version}" | sed -e 's/[^.]*\.\([^.}]*\).*/\1/'`
  if test \! \( \( ${major} -eq ${major_version} \) -a \( ${minor} -ge 6 \) \) ; then
    return 1
  else
    return 0
  fi
}

check_item_compatibility () {
  local item_name="$1"
  local query="query={\"kimcode\":\"${item_name}\"}"
  query="${query}"'&fields={"kim-api-version":1}'
  query="${query}"'&database=obj&history=on'
  local version=`wget -q -O - --post-data="${query}" https://query.openkim.org/api \
                 | \
                 sed -e 's/\[//g' -e 's/\]//g' \
                 -e 's/{"kim-api-version": "\([0-9.]*\)"/\1/g'`
  if test x"" = x"${version}"; then
    printf "*** ERROR *** ${item_name} not found at openkim.org.\n"
    return 1
  else
    if check_version_compatibility "${version}"; then
      return 0
    else
      printf "*** ERROR *** ${item_name} found at openkim.org is not compatible with this version of the KIM API.\n"
      return 1
    fi
  fi
}

check_for_local_build () {
  local name="$1"
  local_build_path=""
  local base_name=`printf -- "${name}" | sed 's|.*/\([^/][^/]*\)/*$|\1|'`
  if test x"${base_name}" != x"${name}"; then
    local_build_path=`printf -- "${name}" | sed 's|/*$||'`
    if test x"/" != x`expr ${local_build_path} : '\(.\)'`; then
      local_build_path="${PWD}/${local_build_path}"
    fi
    return 0
  else
    return 1
  fi
  # output is: local_build_path
}

get_local_build_item_name () {
  local build_path="$1"
  if test \! -d "${build_path}"; then
    printf "Item path '${build_path}' does not exist.\n"
    return 1
  fi

  # create private temporary directory
  if test x"" = x"${TMPDIR}"; then TMPDIR="/tmp"; fi
  local build_dir=`mktemp -d "${TMPDIR}/kim-api-v${major_version}-build-XXXXXXXXXX"`
  if test $? -ne 0; then
    printf "Unable to create temporary directory.\n"
    return 1;
  fi

  local orig_pwd="${PWD}"
  local error="false"
  cd "${build_dir}"

  # setup kim-api config
  ${build_config} --makefile-kim-config > Makefile.KIM_Config

  cp -r "${build_path}" "item"
  cd "item"

  local item_type=`make kim-item-type`
  case "${item_type}" in
    Model|ParameterizedModel|SimulatorModel)
      item_name=`make model-name`
      ;;
    ModelDriver)
      item_name=`make model-driver-name`
      ;;
    *)
      item_name=""
      printf "Item type '${item_type}' unknown.\n"
      error="true"
      ;;
  esac

  cd "${orig_pwd}"
  rm -rf "${build_dir}"
  if test x"true" = x"${error}"; then
    return 1
  else
    return 0
  fi
  # output is: item_name
}

check_makefile_for_v2_build_config () {
  # Assumes we are in the TMP build dir with "${item_id}" as a subdirectory
  local item_id="$1"

  printf "include ./${item_id}/Makefile\n"  > Makefile
  printf 'print-%%: ; @printf -- "$($*)"'  >> Makefile

  local item_build_config=`make print-KIM_API_BUILD_CONFIG 2> /dev/null`

  if test x"kim-api-v2-build-config" = x"${item_build_config}"; then
    printf "Item is designed for kim-api-v2 (Makefile KIM_API_BUILD_CONFIG = kim-api-v2-build-config).\n"
    return 0
  else
    return 1
  fi
}

check_config_file () {
  local config_file_name=`${collections_info} config_file name`
  local drivers_dir=`${collections_info} config_file model_drivers`
  local models_dir=`${collections_info} config_file models`

  if test \! -f "${config_file_name}" -o x"" = x"${drivers_dir}" -o x"" = x"${models_dir}"; then
    printf "Invalid kim-api configuration file.\n"
    return 1
  fi
}

rewrite_config_file_models_dir () {
  if test -d "$1"; then
     local config_file_name=`${collections_info} config_file name`
     local drivers_dir=`${collections_info} config_file model_drivers`
     local models_dir=`cd "$1" && pwd`

     printf "model_drivers_dir = %s\n" "${drivers_dir}" >  "${config_file_name}" || return 1
     printf "models_dir = %s\n" "${models_dir}"         >> "${config_file_name}"
  else
    printf "Directory '%s' does not exist.\n" "$1"
    return 1
  fi
}

rewrite_config_file_drivers_dir () {
  if test -d "$1"; then
    local config_file_name=`${collections_info} config_file name`
    local drivers_dir=`cd "$1" && pwd`
    local models_dir=`${collections_info} config_file models`

    printf "model_drivers_dir = %s\n" "${drivers_dir}" >  "${config_file_name}" || return 1
    printf "models_dir = %s\n" "${models_dir}"         >> "${config_file_name}"
  else
    printf "Directory '%s' does not exist." "$1"
    return 1
  fi
}

get_build_install_item () {
  local install_collection="$1"
  local item_name="$2"
  local use_sudo="$3"
  local PASSWORD="$4"
  local found_item=""
  local item_type=""
  local local_build_path=""

  # check for non-empty item_name
  if test x"" = x"${item_name}"; then
    printf "Empty item id.\n"
    return 1
  fi

  # make changes for installing to CWD collection, if necessary
  if test x"${install_collection}" = x"CWD"; then
    export KIM_API_MODEL_DRIVERS_DIR=${PWD}
    export KIM_API_MODELS_DIR=${PWD}
    install_collection="environment"
  fi

  # check for local build and get item name
  if check_for_local_build "${item_name}"; then  # sets local_build_path
    if ! get_local_build_item_name "${local_build_path}"; then  # sets item_name
      # error message already printed
      return 1
    fi
  fi

  # check for existing item
  if test x"OpenKIM" = x"${item_name}"; then
    found_item=""
    item_type="OpenKIM"
  else
    found_item="`${collections_info} model_drivers find "${item_name}"`"
    found_item="${found_item}""`${collections_info} models find "${item_name}"`"
  fi
  if test x"" != x"${found_item}"; then
    local item_collection=`printf -- "${found_item}" | sed -e 's/ .*//'`
    printf "Item '${item_name}' already installed in collection '${item_collection}'.\n"
    if test x"${item_collection}" = x"${install_collection}"; then
      return 0
    else
      return 1
    fi
  fi

  # create private temporary directory
  if test x"" = x"${TMPDIR}"; then TMPDIR="/tmp"; fi
  local build_dir=`mktemp -d "${TMPDIR}/kim-api-v${major_version}-build-XXXXXXXXXX"`
  if test $? -ne 0; then
    printf "Unable to create temporary directory.\n"
    return 1;
  fi

  (  # subshell
    cd "${build_dir}" || return 1

    # setup kim-api config
    ${build_config} --makefile-kim-config > Makefile.KIM_Config || return 1

    # download item (and possibly its driver)
    if test x"OpenKIM" = x"${item_type}"; then
      local query='query={"type":"mo","kim-api-version":{"$regex":"^'"${major_version}"'\\."}}'
      query="${query}"'&fields={"kimcode":1, "kim-api-version":1}'
      query="${query}"'&database=obj'
      local list=`wget -q -O - --post-data="${query}" https://query.openkim.org/api \
                     | \
                     sed -e 's/\[//g' -e 's/\]//g' \
                     -e 's/{"kim-api-version": "\([0-9.]*\)", "kimcode": "\([^"]*\)"},*/\1:\2/g'`
      for version in ${list}; do \
        if check_version_compatibility "${version}"; then
          get_build_install_item "$install_collection" "${modname}" "${use_sudo}" "${PASSWORD}" || return 1
        fi
      done
    else
      if test x"" = x"${local_build_path}"; then
        if check_item_compatibility "${item_name}"; then
          printf "Downloading..............@%s\n" "${item_name}" | sed -e 's/ /./g' -e 's/@/ /g'
          if wget -q --content-disposition "https://openkim.org/download/${item_name}.txz"; then
            tar Jxf "${item_name}.txz" 2>&1 | sed -e 's/^/                /' &&
              rm -f "${item_name}.txz"
          else
            printf "                Unable to download ${item_name} from https://openkim.org.  Check the KIM Item ID for errors.\n"
            return 1
          fi
        else
          return 1
        fi
      else
        cp -r "${local_build_path}" "./${item_name}"
      fi
      cd "./${item_name}"
      item_type="`${make_command} kim-item-type`"
      if (cd .. && check_makefile_for_v2_build_config "${item_name}"); then
        return 1
      elif test 0 -lt `grep -c MAKE_SYSTEM Makefile`; then
        printf "*** ERROR *** ${item_name} appears to be written for an older, incompatible, version of the KIM API.\n"
        return 1
      elif test x"ParameterizedModel" = x"${item_type}"; then
        dvr="`${make_command} model-driver-name`"
        if test x"" != x"`${collections_info} model_drivers find "${dvr}"`"; then
          printf "Using@installed@driver...@%s\n" "${dvr}" | sed -e 's/ /./g' -e 's/@/ /g' || return 1
        else
          printf "This model needs its driver '${dvr}', trying to find it at openkim.org.\n"
          # try openkim.org first
          if ! get_build_install_item "${install_collection}" "${dvr}" "${use_sudo}" "${PASSWORD}"; then
            if test x"" != x"${local_build}"; then
              # now try local
              printf "Now trying to find '${dvr}' locally.\n"
              dvr=`printf "${local_build}" | sed "s|^\(.*/\)${item_name}\$|\1${dvr}|"`
              get_build_install_item "${install_collection}" "${dvr}" "${use_sudo}" "${PASSWORD}" || return 1
            fi
          fi
        fi
        if ! ((${make_command} clean && ${make_command}) > ./make-log.txt 2>&1 && \
                if test x"sudo-yes" = x"${use_sudo}"; then
                  printf -- "${PASSWORD}\n" | sudo -k -S ${make_command} "install-${install_collection}" 2> /dev/null
                else
                  ${make_command} "install-${install_collection}"
                fi); then
          cat ./make-log.txt
          return 1
        fi
      elif test x"Model"          = x"${item_type}" -o \
                x"ModelDriver"    = x"${item_type}" -o \
                x"SimulatorModel" = x"${item_type}"; then
        if ! ((${make_command} clean && ${make_command}) > ./make-log.txt 2>&1 && \
                if test x"sudo-yes" = x"${use_sudo}"; then
                  printf -- "${PASSWORD}\n" | sudo -k -S ${make_command} "install-${install_collection}" 2> /dev/null
                else
                  ${make_command} "install-${install_collection}"
                fi); then
          cat ./make-log.txt
          return 1
        fi
      else
        printf "Item '${item_name}' of unknown type '${item_type}'.\n"
        return 1
      fi
    fi
  ) || return 1  # exit subshell

  rm -rf "${build_dir}" || return 1
}

remove_item () {
  local item_name="$1"
  local use_sudo="$2"
  local PASSWORD="$3"
  local found_item=""
  local item_type=""

  # check for existing item
  found_item="`${collections_info} model_drivers find "${item_name}"`"
  if test x"" = x"${found_item}"; then
    found_item="`${collections_info} models find "${item_name}"`"
    if test x"" = x"${found_item}"; then
      printf "Item '${item_name}' not installed.\n"
      return 1
    else
      item_type="models"
    fi
  else
    item_type="model_drivers"
  fi

  local item_dir=`${collections_info} "${item_type}" find "${item_name}" | sed -e 's/^[^ ]* [^ ]* \([^ ]*\).*/\1/'`"/${item_name}"

  printf "Removing '%s'.\n" "${item_dir}"
  if test x"sudo-yes" = x"${use_sudo}"; then
    printf -- "${PASSWORD}\n" | sudo -k -S rm -rf "${item_dir}" 2> /dev/null || return 1
  else
    rm -rf "${item_dir}" || return 1
  fi
}

split_drivers_list_into_collections () {
  local with_version=$2

  drivers_cwd_collection=""; number_drivers_cwd=0
  drivers_env_collection=""; number_drivers_env=0
  drivers_usr_collection=""; number_drivers_usr=0
  drivers_sys_collection=""; number_drivers_sys=0
    while read line; do
      local collection=`printf -- "$line" | sed -e 's/\([^ ]*\) .*/\1/'`
      local name=`printf -- "$line" | sed -e 's/[^ ]* \([^ ]*\) .*/\1/'`
      local directory=`printf -- "$line" | sed -e 's/[^ ]* [^ ]* \([^ ]*\) .*/\1/'`
      local version=`printf -- "$line" | sed -e 's/[^ ]* [^ ]* [^ ]* \([^ ]*\).*/\1/'`
      case $collection in
        "")
        # empty do nothing
        ;;
        CWD)
          number_drivers_cwd=`expr $number_drivers_cwd \+ 1`
          drivers_cwd_collection="${drivers_cwd_collection}\t${name}"
          if test x"${with_version}" = x"yes"; then
            drivers_cwd_collection="${drivers_cwd_collection}\n\t\t${version}"
          fi
          drivers_cwd_collection="${drivers_cwd_collection}\n"
          ;;
        environment)
          number_drivers_env=`expr $number_drivers_env \+ 1`
          drivers_env_collection="${drivers_env_collection}\t${name}"
          if test x"${with_version}" = x"yes"; then
            drivers_env_collection="${drivers_env_collection}\n\t\t${version}"
          fi
          drivers_env_collection="${drivers_env_collection}\n"
          ;;
        user)
          number_drivers_usr=`expr $number_drivers_usr \+ 1`
          drivers_usr_collection="${drivers_usr_collection}\t${name}"
          if test x"${with_version}" = x"yes"; then
            drivers_usr_collection="${drivers_usr_collection}\n\t\t${version}"
          fi
          drivers_usr_collection="${drivers_usr_collection}\n"
          ;;
        system)
          number_drivers_sys=`expr $number_drivers_sys \+ 1`
          drivers_sys_collection="${drivers_sys_collection}\t${name}"
          if test x"${with_version}" = x"yes"; then
            drivers_sys_collection="${drivers_sys_collection}\n\t\t${version}"
          fi
          drivers_sys_collection="${drivers_sys_collection}\n"
          ;;
        *)
          printf "Error unknown collection!\n"
          exit 1
          ;;
      esac
    done <<EOF
$1
EOF
}

split_models_list_into_collections () {
  local with_version=$2

  models_cwd_collection=""; number_models_cwd=0
  models_env_collection=""; number_models_env=0
  models_usr_collection=""; number_models_usr=0
  models_sys_collection=""; number_models_sys=0
    while read line; do
      local collection=`printf -- "$line" | sed -e 's/\([^ ]*\) .*/\1/'`
      local name=`printf -- "$line" | sed -e 's/[^ ]* \([^ ]*\) .*/\1/'`
      local directory=`printf -- "$line" | sed -e 's/[^ ]* [^ ]* \([^ ]*\) .*/\1/'`
      local version=`printf -- "$line" | sed -e 's/[^ ]* [^ ]* [^ ]* \([^ ]*\).*/\1/'`
      case $collection in
        "")
        # empty do nothing
        ;;
        CWD)
          number_models_cwd=`expr $number_models_cwd \+ 1`
          models_cwd_collection="${models_cwd_collection}\t${name}"
          if test x"${with_version}" = x"yes"; then
            models_cwd_collection="${models_cwd_collection}\n\t\t${version}"
          fi
          models_cwd_collection="${models_cwd_collection}\n"
          ;;
        environment)
          number_models_env=`expr $number_models_env \+ 1`
          models_env_collection="${models_env_collection}\t${name}"
          if test x"${with_version}" = x"yes"; then
            models_env_collection="${models_env_collection}\n\t\t${version}"
          fi
          models_env_collection="${models_env_collection}\n"
          ;;
        user)
          number_models_usr=`expr $number_models_usr \+ 1`
          models_usr_collection="${models_usr_collection}\t${name}"
          if test x"${with_version}" = x"yes"; then
            models_usr_collection="${models_usr_collection}\n\t\t${version}"
          fi
          models_usr_collection="${models_usr_collection}\n"
          ;;
        system)
          number_models_sys=`expr $number_models_sys \+ 1`
          models_sys_collection="${models_sys_collection}\t${name}"
          if test x"${with_version}" = x"yes"; then
            models_sys_collection="${models_sys_collection}\n\t\t${version}"
          fi
          models_sys_collection="${models_sys_collection}\n"
          ;;
        *)
          printf "Error unknown collection!\n"
          exit 1
          ;;
      esac
    done <<EOF
$1
EOF
}

print_separator_line () {
  printf "%79s\n" " " | sed -e "s/ /$1/g"
}

print_list_of_drivers_and_models () {
  local with_version=$1

  config_env_name=`${collections_info} config_file env | sed -e 's/ .*//'`
  config_env=`${collections_info} config_file env | sed -e 's/^[^ ]* //'`
  if test x"" = x"${config_env}"; then config_env="--empty--"; fi
  drivers_env_name=`${collections_info} env env | sed -e 's/^[^ ]* //'`
  drivers_env=`${collections_info} env model_drivers`
  if test x"" = x"${drivers_env}"; then drivers_env="--empty--"; fi
  models_env_name=`${collections_info} env env | sed -e 's/ .*//'`
  models_env=`${collections_info} env models`
  if test x"" = x"${models_env}"; then models_env="--empty--"; fi

  printf "\n\n"
  printf "Knowledgebase of Interatomic Models (KIM)"
  printf -- "  ---  Model Collections Listing\n"
  print_separator_line "="
  printf "\n"
  printf "kim-api library: \n\t%s\n" `${collections_info} system library`
  printf "\n"
  printf "kim-api configuration file:\n\t%s\n" \
         `${collections_info} config_file name`
  printf "\n\n"
  printf "Environment Variables:\n"
  print_separator_line "-"
  printf -- "${config_env_name}:\n"
  printf -- "\t${config_env}\n"
  printf "\n"
  printf -- "${drivers_env_name}:\n"
  printf -- "%s\n" "`printf -- "${drivers_env}" | sed -e 's/^/	/g'`"
  printf "\n"
  printf -- "${models_env_name}:\n"
  printf -- "%s\n" "`printf -- "${models_env}" | sed -e 's/^/	/g'`"
  printf "\n"
  print_separator_line "="


  model_drivers_list=`${collections_info} model_drivers`
  split_drivers_list_into_collections "${model_drivers_list}" "${with_version}"
  models_list=`${collections_info} models`
  split_models_list_into_collections "${models_list}" "${with_version}"

  printf "\n\n\n"
  printf "Current Working Directory Collection\n"
  print_separator_line "-"
  printf "Drivers: %s\n" "${PWD}"
  if test $number_drivers_cwd -gt 0; then
    printf "${drivers_cwd_collection}"
  else
    printf "\t--empty--\n"
  fi
  printf "\n"

    printf "Models: %s\n" "${PWD}"
  if test $number_models_cwd -gt 0; then
    printf -- "${models_cwd_collection}"
  else
    printf "\t--empty--\n"
  fi
  printf "\n\n"

  printf "Environment Variable Collection\n"
  print_separator_line "-"
  printf "Drivers: "
  if test x"--empty--" = x"${drivers_env}"; then
    printf "%s" ${drivers_env}
  else
    printf "'%s' " ${drivers_env}
  fi
  printf "\n"
  if test $number_drivers_env -gt 0; then
    printf -- "${drivers_env_collection}"
  else
    printf "\t--empty--\n"
  fi
  printf "\n"
  printf "Models: "
  if test x"--empty--" = x"${models_env}"; then
    printf "%s" ${models_env}
  else
    printf "'%s' " ${models_env}
  fi
  printf "\n"
  if test $number_models_env -gt 0; then
    printf -- "${models_env_collection}"
  else
    printf "\t--empty--\n"
  fi
  printf "\n\n"

  drivers_usr=`${collections_info} config_file model_drivers`
  if test x"" = x"${drivers_usr}"; then drivers_usr="--empty--"; fi
  models_usr=`${collections_info} config_file models`
  if test x"" = x"${models_usr}"; then models_usr="--empty--"; fi
  printf "User Collection\n"
  print_separator_line "-"
  printf "Drivers: %s\n" "${drivers_usr}"
  if test $number_drivers_usr -gt 0; then
    printf -- "${drivers_usr_collection}"
  else
    printf "\t--empty--\n"
  fi
  printf "\n"
  printf "Models: %s\n" "${models_usr}"
  if test $number_models_usr -gt 0; then
    printf -- "${models_usr_collection}"
  else
    printf "\t--empty--\n"
  fi
  printf "\n\n"

  printf "System Collection\n"
  print_separator_line "-"
  printf "Drivers: %s\n" `${collections_info} system model_drivers`
  if test $number_drivers_sys -gt 0; then
    printf -- "${drivers_sys_collection}"
  else
    printf "\t--empty--\n"
  fi
  printf "\n\n"
  printf "Models: %s\n" `${collections_info} system models`
  if test $number_models_sys -gt 0; then
    printf -- "${models_sys_collection}"
  else
    printf "\t--empty--\n"
  fi
  printf "\n"
}


get_password () {
  printf "Enter Password : "
  stty -echo
  trap 'stty echo' EXIT
  read PASSWORD
  stty echo
  trap - EXIT
  printf "\n"
  if ! (printf -- "${PASSWORD}\n" | \
          sudo -k -S printf "" > /dev/null 2>&1); then
    printf "Bad password.\n"
    return 1
  fi
}

get_confirmation () {
  local ANSWER
  printf "[y/n] : "
  read ANSWER
  if test x"${ANSWER}" = x"y"; then
    return 0;
  else
    return 1;
  fi
}

######## main script ########

# check that command is given
if test $# -lt 1; then
  usage
  exit 1
else
  command=$1
  case $command in
    list|set-user-drivers-dir|set-user-models-dir|install|reinstall|remove|remove-all)
    ;;
    *)
      printf "unknown command: %s\n\n" $command
      usage
      exit 1
  esac
fi

if ! check_config_file; then
  printf "Aborting!\n"
  exit 1
fi

case $command in
  list)
    if test $# -gt 2; then
      usage
      exit 1
    elif test x"$2" = x"--with-version"; then
      print_list_of_drivers_and_models "yes"
    else
      print_list_of_drivers_and_models "no"
    fi
    ;;
  set-user-models-dir)
    if test $# -lt 2; then
      usage
      exit 1
    else
      subcommand=$2
      if ! rewrite_config_file_models_dir "$subcommand"; then
        printf "\nAborting!\n"
        exit 1
      else
        printf "\nSuccess!\n."
      fi
    fi
    ;;
  set-user-drivers-dir)
    if test $# -lt 2; then
      usage
      exit 1
    else
      subcommand=$2
      if ! rewrite_config_file_drivers_dir "$subcommand"; then
        printf "\nAborting!\n"
        exit 1
      else
        printf "\nSuccess!\n"
      fi
    fi
    ;;
  install)
    if test $# -lt 3; then
      usage
      exit 1
    else
      subcommand=$2
      shift
      shift
      case $subcommand in
        CWD|environment|user)
          for item_id in $@; do
            if ! get_build_install_item "${subcommand}" "${item_id}" "sudo-no" ""; then
              printf "\nAborting!\n"
              exit 1
            else
              printf "\nSuccess!\n"
            fi
          done
          ;;
        system)
          PASSWORD=""
          if test x"--sudo" = x"$1"; then
            shift
            use_sudo="sudo-yes"
            if ! get_password; then
              printf "\nAborting!\n"
              exit 1
            fi
          else
            use_sudo="sudo-no"
          fi
          for item_id in $@; do
            if ! get_build_install_item "system" "${item_id}" "${use_sudo}" "${PASSWORD}"; then
              printf "\nAborting!\n"
              exit 1
            else
              printf "\nSuccess!\n"
            fi
          done
          ;;
        *)
          printf "unknown subcommand: %s\n\n" $subcommand
          usage
          exit 1
          ;;
      esac
    fi
    ;;
  reinstall)
    if test $# -lt 2; then
      usage
      exit 1
    else
      shift
      if test x"--sudo" = x"$1"; then
        shift
        use_sudo="sudo-yes"
        if ! get_password; then
          printf "\nAborting!\n"
          exit 1
        fi
      else
        use_sudo="sudo-no"
      fi
    fi
    for item_id in $@; do
      if check_for_local_build "${item_id}"; then  # sets local_build_path
        if ! get_local_build_item_name "${local_build_path}"; then  # sets item_name
          printf "\nAborting!\n"
          exit 1
        fi
      else
        if ! check_item_compatibility "${item_id}"; then
          printf "\nAborting!\n"
          exit 1
        else
          item_name="${item_id}"
        fi
      fi
      found_item="`${collections_info} model_drivers find "${item_name}"` `${collections_info} models find "${item_name}"`"
      if test \! x"" = x"${found_item}"; then
        item_collection=`printf "${found_item}" | sed -e 's/^[[:space:]]*\([^[:space:]]*\) .*/\1/'`
        if ! (remove_item "${item_name}" "${use_sudo}" "${PASSWORD}" && \
                get_build_install_item "${item_collection}" "${item_id}" "${use_sudo}" "${PASSWORD}"); then
          printf "\nAborting!\n"
          exit 1
        else
          printf "\nSuccess!\n"
        fi
      else
        printf "\nAborting!\n"
        exit 1
      fi
    done
    ;;
  remove)
    if test $# -lt 2; then
      usage
      exit 1
    else
      shift
      if test x"--sudo" = x"$1"; then
        shift
        use_sudo="sudo-yes"
        if ! get_password; then
          printf "\nAborting!\n"
          exit 1
        fi
      else
        use_sudo="sudo-no"
      fi
    fi
    printf "This will remove all files associated with these items.\n"
    printf "\n"
    printf "Are you sure you want to proceed? "
    if get_confirmation; then
      for item_id in $@; do
        if ! remove_item "${item_id}" "${use_sudo}" "${PASSWORD}"; then
          printf "\nAborting!\n"
          exit 1
        else
          printf "\nSuccess!\n"
        fi
      done
    else
      printf "\nAborting!\n"
    fi
    ;;
  remove-all)
    if test $# -eq 2; then
      if test x"--sudo" = x"$2"; then
        use_sudo="sudo-yes"
        if ! get_password; then
          printf "\nAborting!\n"
          exit 1
        fi
      else
        use_sudo="sudo-no"
      fi
    fi
    printf "This will remove all files associated with all items in the 'environment',\n"
    printf "'user', and 'system' collections.\n"
    printf "\n"
    printf "Are you sure you want to proceed? "
    if get_confirmation; then
      # get full list and loop over it to remove
      for item_id in `${collections_info} model_drivers | sed -e 's/^[^[:space:]]* \([^[:space:]]*\).*/\1/'`; do
        if ! remove_item "${item_id}" "${use_sudo}" "${PASSWORD}"; then
          printf "\n Aborting!\n"
          exit 1
        fi
      done
      for item_id in `${collections_info} models | sed -e 's/^[^[:space:]]* \([^[:space:]]*\).*/\1/'`; do
        if ! remove_item "${item_id}" "${use_sudo}" "${PASSWORD}"; then
          printf "\n Aborting!\n"
          exit 1
        fi
      done
      printf "\nSuccess!\n"
    else
      printf "\nAborting!\n"
      exit 1
    fi
    ;;
esac
