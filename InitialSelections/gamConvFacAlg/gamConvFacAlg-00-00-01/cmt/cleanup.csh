# echo "cleanup gamConvFacAlg gamConvFacAlg-00-00-01 in /afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtgamConvFacAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtgamConvFacAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=gamConvFacAlg -version=gamConvFacAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  $* >${cmtgamConvFacAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=gamConvFacAlg -version=gamConvFacAlg-00-00-01 -path=/afs/ihep.ac.cn/users/y/yuanxq/higgs/workArea/704  $* >${cmtgamConvFacAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtgamConvFacAlgtempfile}
  unset cmtgamConvFacAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtgamConvFacAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtgamConvFacAlgtempfile}
unset cmtgamConvFacAlgtempfile
exit $cmtcleanupstatus

