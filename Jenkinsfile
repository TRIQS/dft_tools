def projectName = "app4triqs" /* set to app/repo name */

/* which platform to build documentation on */
def documentationPlatform = "ubuntu-clang"
/* depend on triqs upstream branch/project */
def triqsBranch = env.CHANGE_TARGET ?: env.BRANCH_NAME
def triqsProject = '/TRIQS/triqs/' + triqsBranch.replaceAll('/', '%2F')
/* whether to publish the results (disabled for template project) */
def publish = !env.BRANCH_NAME.startsWith("PR-") && projectName != "app4triqs"

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

/* map of all builds to run, populated below */
def platforms = [:]

/****************** linux builds (in docker) */
/* Each platform must have a cooresponding Dockerfile.PLATFORM in triqs/packaging */
def dockerPlatforms = ["ubuntu-clang", "ubuntu-gcc", "centos-gcc"]
/* .each is currently broken in jenkins */
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
      def img = docker.build("flatironinstitute/${projectName}:${env.BRANCH_NAME}-${env.STAGE_NAME}", "--build-arg BUILD_DOC=${platform==documentationPlatform} .")
      if (!publish || platform != documentationPlatform) {
        /* but we don't need the tag so clean it up (except for documentation) */
        sh "docker rmi --no-prune ${img.imageName()}"
    }
  } }
  } }
}

/****************** osx builds (on host) */
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
      def triqsDir = "${env.HOME}/install/triqs/${triqsBranch}/${platform}"
      dir(installDir) {
        deleteDir()
      }

      checkout scm
      dir(buildDir) { withEnv(platformEnv[1]+[
        "PATH=$triqsDir/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin",
        "CPATH=$triqsDir/include",
        "LIBRARY_PATH=$triqsDir/lib",
        "CMAKE_PREFIX_PATH=$triqsDir/share/cmake"]) {
        deleteDir()
      sh "cmake $srcDir -DCMAKE_INSTALL_PREFIX=$installDir -DTRIQS_ROOT=$triqsDir"
      sh "make -j3"
        try {
          sh "make test"
        } catch (exc) {
          archiveArtifacts(artifacts: 'Testing/Temporary/LastTest.log')
          throw exc
        }
        sh "make install"
      } }
    } }
  } }
}

/****************** wrap-up */
try {
  parallel platforms
  if (publish) { node("docker") {
    /* Publish results */
    stage("publish") { timeout(time: 1, unit: 'HOURS') {
      def commit = sh(returnStdout: true, script: "git rev-parse HEAD").trim()
      def workDir = pwd()
      /* Update documention on gh-pages branch */
      dir("$workDir/gh-pages") {
        def subdir = env.BRANCH_NAME
        git(url: "ssh://git@github.com/TRIQS/${projectName}.git", branch: "gh-pages", credentialsId: "ssh", changelog: false)
        sh "rm -rf ${subdir}"
        docker.image("flatironinstitute/${projectName}:${env.BRANCH_NAME}-${documentationPlatform}").inside() {
          sh "cp -rp \$INSTALL/share/doc/${projectName} ${subdir}"
        }
        sh "git add -A ${subdir}"
        sh """
          git commit --author='Flatiron Jenkins <jenkins@flatironinstitute.org>' --allow-empty -m 'Generated documentation for ${env.BRANCH_NAME}' -m '${env.BUILD_TAG} ${commit}'
        """
        // note: credentials used above don't work (need JENKINS-28335)
        sh "git push origin gh-pages"
      }
      /* Update docker repo submodule */
      dir("$workDir/docker") { try {
        git(url: "ssh://git@github.com/TRIQS/docker.git", branch: env.BRANCH_NAME, credentialsId: "ssh", changelog: false)
        sh "echo '160000 commit ${commit}\t${projectName}' | git update-index --index-info"
        sh """
          git commit --author='Flatiron Jenkins <jenkins@flatironinstitute.org>' -m 'Autoupdate ${projectName}' -m '${env.BUILD_TAG}'
        """
        // note: credentials used above don't work (need JENKINS-28335)
        sh "git push origin ${env.BRANCH_NAME}"
      } catch (err) {
	/* Ignore, non-critical -- might not exist on this branch */
        echo "Failed to update docker repo"
      } }
    } }
  } }
} catch (err) {
  /* send email on build failure (declarative pipeline's post section would work better) */
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
    to: 'nils.wentzell@gmail.com, dsimon@flatironinstitute.org',
    recipientProviders: [
      [$class: 'DevelopersRecipientProvider'],
    ],
    replyTo: '$DEFAULT_REPLYTO'
  )
  throw err
}
