def projectName = "dft_tools"
def triqsBranch = env.CHANGE_TARGET ?: env.BRANCH_NAME
def triqsProject = '/TRIQS/triqs/' + triqsBranch.replaceAll('/', '%2F')
def documentationPlatform = "ubuntu-clang"

properties([
  disableConcurrentBuilds(),
  buildDiscarder(logRotator(numToKeepStr: '10', daysToKeepStr: '30')),
  pipelineTriggers([
    upstream(
      threshold: 'SUCCESS',
      upstreamProjects: triqsProject
    )
  ])
])

def platforms = [:]

def dockerPlatforms = ["ubuntu-clang", "ubuntu-gcc", "centos-gcc"]
for (int i = 0; i < dockerPlatforms.size(); i++) {
  def platform = dockerPlatforms[i]
  platforms[platform] = { -> node('docker') {
    stage(platform) { timeout(time: 1, unit: 'HOURS') {
      checkout scm
      /* construct a Dockerfile for this base */
      sh """
        ( echo "FROM flatironinstitute/triqs:${triqsBranch}-${env.STAGE_NAME}" ; sed '0,/^FROM /d' Dockerfile ) > Dockerfile.jenkins
        mv -f Dockerfile.jenkins Dockerfile
      """
      /* build and tag */
      def img = docker.build("flatironinstitute/${projectName}:${env.BRANCH_NAME}-${env.STAGE_NAME}", "--build-arg Build_Documentation=${platform==documentationPlatform} .")
      if (env.BRANCH_NAME.startsWith("PR-") || platform != documentationPlatform) {
        /* but we don't need the tag so clean it up (except for documentation) */
        sh "docker rmi ${img.imageName()}"
      }
    } }
  } }
}

def osxPlatforms = [
  ["gcc", ['CC=gcc-7', 'CXX=g++-7']],
  ["clang", ['CC=/usr/local/opt/llvm/bin/clang', 'CXX=/usr/local/opt/llvm/bin/clang++', 'CXXFLAGS=-I/usr/local/opt/llvm/include', 'LDFLAGS=-L/usr/local/opt/llvm/lib']]
]
for (int i = 0; i < osxPlatforms.size(); i++) {
  def platformEnv = osxPlatforms[i]
  def platform = platformEnv[0]
  platforms["osx-$platform"] = { -> node('osx && triqs') {
    stage("osx-$platform") { timeout(time: 1, unit: 'HOURS') {
      def srcDir = pwd()
      def tmpDir = pwd(tmp:true)
      def buildDir = "$tmpDir/build"
      def installDir = "$tmpDir/install"

      dir(installDir) {
        deleteDir()
      }

      copyArtifacts(projectName: triqsProject, selector: upstream(fallbackToLastSuccessful: true), filter: "osx-${platform}.zip")
      unzip(zipFile: "osx-${platform}.zip", dir: installDir)
      /* fixup zip-stripped permissions (JENKINS-13128) */
      sh "chmod +x $installDir/bin/*"

      checkout scm

      dir(buildDir) { withEnv(platformEnv[1]+[
          "PATH=$installDir/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin",
          "CPATH=$installDir/include",
          "LIBRARY_PATH=$installDir/lib",
          "CMAKE_PREFIX_PATH=$installDir/share/cmake"]) {
        deleteDir()
        sh "cmake $srcDir -DTRIQS_ROOT=$installDir"
        sh "make -j2"
        try {
          sh "make test"
        } catch (exc) {
          archiveArtifacts(artifacts: 'Testing/Temporary/LastTest.log')
          throw exc
        }
        sh "make install"
      } }
      // zip(zipFile: "osx-${platform}.zip", archive: true, dir: installDir)
    } }
  } }
}

try {
  parallel platforms
  if (!env.BRANCH_NAME.startsWith("PR-")) {
    node("docker") {
      stage("documentation") { timeout(time: 1, unit: 'HOURS') {
        def workDir = pwd()
        dir("$workDir/gh-pages") {
          def subdir = env.BRANCH_NAME
          git(url: "ssh://git@github.com/TRIQS/${projectName}.git", branch: "gh-pages", credentialsId: "ssh", changelog: false)
          sh "rm -rf ${subdir}"
          docker.image("flatironinstitute/${projectName}:${env.BRANCH_NAME}-${documentationPlatform}").inside() {
            sh "cp -rp \$INSTALL/share/doc/${projectName} ${subdir}"
          }
          sh "git add -A ${subdir}"
          sh """
            git commit --author='Flatiron Jenkins <jenkins@flatironinstitute.org>' --allow-empty -m 'Generated documentation for ${env.BRANCH_NAME}' -m "`git --git-dir ${workDir}/.git rev-parse HEAD`"
          """
          // note: credentials used above don't work (need JENKINS-28335)
          sh "git push origin gh-pages"
        }
      } }
    }
  }
} catch (err) {
  if (env.BRANCH_NAME != "jenkins") emailext(
    subject: "\$PROJECT_NAME - Build # \$BUILD_NUMBER - FAILED",
    body: """\$PROJECT_NAME - Build # \$BUILD_NUMBER - FAILED

$err

Check console output at \$BUILD_URL to view full results.

Building \$BRANCH_NAME for \$CAUSE
\$JOB_DESCRIPTION

Chages:
\$CHANGES

End of build log:
\${BUILD_LOG,maxLines=60}
    """,
    to: 'mzingl@flatironinstitute.org, hstrand@flatironinstitute.org, nils.wentzell@gmail.com, dsimon@flatironinstitute.org',
    recipientProviders: [
      [$class: 'DevelopersRecipientProvider'],
    ],
    replyTo: '$DEFAULT_REPLYTO'
  )
  throw err
}
