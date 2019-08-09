# echo "setup gamConvFacAlg gamConvFacAlg-00-00-01 in /afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtgamConvFacAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtgamConvFacAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=gamConvFacAlg -version=gamConvFacAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  -no_cleanup $* >${cmtgamConvFacAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=gamConvFacAlg -version=gamConvFacAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  -no_cleanup $* >${cmtgamConvFacAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtgamConvFacAlgtempfile}
  unset cmtgamConvFacAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtgamConvFacAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtgamConvFacAlgtempfile}
unset cmtgamConvFacAlgtempfile
exit $cmtsetupstatus

